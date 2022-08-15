# Dependencies 

library(data.table)         
library(tidyverse)
library(Rcpp)
library(parallel)
library(flexmix)
library(mhsmm)
library(umap)
library(round)

setDTthreads(threads = 0, restore_after_fork = T)

sourceCpp("/hpc/hub_oudenaarden/vincentvb/forkdot/src/forkdot_src.cpp")

mstep.exp <- function(x, wt)       {list(lambda = apply(wt, 2, function(w){weighted.mean(x, w)}))}
rexp.hsmm <- function(j, model)    {rexp(1, 1 / model$parms.emission$lambda[j])}
dexp.hsmm <- function(x, j, model) {dexp(x, 1 / model$parms.emission$lambda[j])}

# processing functions

super_poisson_test.files <- function(files){
 
  cuts_cell <- rbindlist(map(files, function(x){
    cat("=")
    fread(  x) %>%
      {if(ncol(.) == 5){.[,V6 := "none"][]}else{.}} %>%
      setNames(   c("cell", "chr", "bp", "strand", "readlength", "allele")) %>%
      filter(   chr %in% c(1:22, "X"))
  }))[, bp_round := as.integer(round((bp + 49999.99) / 100000) * 100000 - 49999.99)
  ][,.(N = as.double(.N)), .(cell, chr, bp_round)
  ][,.(mean = sum(N) / 29111,
       variance = sum((N - (sum(N) / 29111))^2) / 29111), .(cell)
  ][, poisson_score := log(sqrt(variance) / mean) + 0.5 * log(mean)]
  
  return(cuts_cell)
}

super_filterer.files     <- function(experiment, files){
  
  print(experiment)
  
  cuts_cell <- rbindlist(map(files, function(x){
    cat("=")
    plate <- fread(x)
    if(ncol(plate) == 6){
    }else if(is.numeric(plate$V5)){
      plate[, V6 := "none"][]} else {plate[,V6 := V5];
        plate[, V5 := 42L][]}
    
    setNames(plate, c("cell", "chr", "bp", "strand", "readlength", "allele"))[chr %in% c(1:22, "X")]
  }))
  
  cuts_cell[,readlength := as.numeric(readlength)]
  
  #head(cuts_cell, 10)
  
  poisson_test <- fread(paste0("hpc/vincentvb/forkdot/poisson_test/", experiment, ".tsv")
  )[order(cell),cell2 := 1:.N
  ][poisson_score > 0.1 
  ][log(mean) > ifelse(str_detect(experiment, "JvB-exp049-15m-EdU-beads"), -1, -2.5)
  ][log(mean) < ifelse(str_detect(experiment, "JvB-exp049-15m-EdU-beads"), 1, 2.5)]
  
  exp_fit_cell.function <- function(single_cell){
    
    sc  <- single_cell[,dist := as.integer(c(0, diff(bp))), .(chr)]
    
    YYY <- flexmix(dist ~ 1, model = list(FLXMCdist1(dist = "exp")), k = 2, 
                   data = sc[dist > 0])
    
    merge.data.table(sc[,.(chr, bp, allele)], 
                     cbind(sc[dist > 0][,.(bp, chr)], 
                           YYY@posterior$scaled[,which.max(parameters(YYY))]),
                     all.x = T, by = c("bp", "chr")
    )[order(bp), posterior := pmax(V2, shift(V2, -1), na.rm = T), .(chr)
    ][, V2 := NULL]
    
  }
  
  exp_fit_cell <- cuts_cell[(cell %in% poisson_test$cell) & (readlength < 1000)
  ][, cell := setNames(poisson_test$cell2, poisson_test$cell)[cell]
  ] %>% split(.$cell) %>%
    mclapply(exp_fit_cell.function, mc.cores = 10) %>%
    rbindlist(idcol = "cell")

  exp_fit_cell[,posterior := round(posterior * 1000)]
  
  return(exp_fit_cell)
}

super_sphase_order       <- function(inputfile){
  
  dist_matrix <- inputfile %>% 
    pivot_wider(names_from = from, values_from = dist) %>%
    column_to_rownames("cell")
  
  umap.custom                <- umap.defaults
  umap.custom$n_neighbors    =  round(log2(nrow(dist_matrix)))
  umap.custom$n_components   =  2
  umap.custom$n_epochs       =  200
  umap.custom$input          =  "dist"
  umap.custom$verbose        =  F
  
  umap.custom2               <- umap.custom
  umap.custom2$n_neighbors   =  round(log2(nrow(dist_matrix)) * 2)
  umap.custom2$n_components  =  1
  
  umap_matrix <- mclapply(1:100, function(x){
    
    set.seed(1000 + x)
    cell_umap2      <- umap(as.matrix(as.dist(dist_matrix)), config = umap.custom2)
    
    scale(cell_umap2$layout)
    
  }, mc.cores = 5)
  
  bootstrap_threshold <- 0.85
  cluster_threshold   <- 0.8
  
  run_select <- colMeans(cor(bind_cols(umap_matrix), method = "spearman")) > bootstrap_threshold
  
  if(sum(run_select) == 0) { print("better check your bootstrap results")
    return( colMeans(cor(bind_cols(umap_matrix), method = "spearman")))}
  
  umap_tour_bootstrap <- as.data.table(
    data.frame(bind_cols(umap_matrix[run_select]), 
               cell = names(dist_matrix)) %>% 
      pivot_longer(-cell)
  )[order(value), cluster := .(cumsum(c(0, diff(value) > 0.1))), .(cell)
  ][, N := .(.N / sum(run_select)), .(cell, cluster)
  ][, N_threshold := .(max(N) > cluster_threshold), .(cell)
  ][, mean := .(weighted.mean(value, N > ifelse(N_threshold[1], cluster_threshold, 0))), .(cell)
  ][, stdev := .(sd(value[N > cluster_threshold])), .(cell)
  ][, rank := .(as.double(dense_rank(mean))), .(N_threshold)
  ][, rank := .(rank / max(rank)), .(N_threshold)] 
  
  selected_cells <- umap_tour_bootstrap[(N_threshold) & !is.na(stdev)] %>% 
    distinct(cell, .keep_all = T)
  
  cell_umap  <- umap(as.matrix(as.dist(dist_matrix[selected_cells$cell, selected_cells$cell])), 
                     config = umap.custom)
  
  umap_tour_final <- data.frame(cell_umap$knn$indexes,
                                cell_umap$layout) %>% 
    rownames_to_column(   "cell") %>% 
    left_join(selected_cells[,.(cell, mean, stdev, rank, N)], by = "cell") %>%
    pivot_longer(c(-cell, -X1.1, -X2.1, -mean, -stdev, -rank, -N)) %>%
    mutate(X2 = cell_umap$layout[value,1],
           Y2 = cell_umap$layout[value,2])
  
  return(list(bootstrap = umap_tour_bootstrap, 
               final_order = umap_tour_final))
}

super_overlap_score      <- function(inputfile){
  
  smooth_chr_cell <- inputfile[posterior > 500, .(cell = as.character(cell), chr, bp)] %>% 
    data.table:::split.data.table(by = "cell", keep.by = F) %>%
    map(function(celler){
      try({
        celler[, round_bp := .(list(as.integer(((round(bp / 5000) - 5):(round(bp / 5000) + 5)) * 5000))), .(chr, bp)
        ][, density  := .(list(dnorm(round_bp[[1]], bp, sd = 25000 / 3))), .(chr, bp)
        ][, .(round_bp = unlist(round_bp), density = unlist(density)), .(chr)
        ][, .(density = sum(density)), .(chr, round_bp)
        ][, tot_density := sum(density), .(chr)] })
    }) %>% rbindlist(idcol = "cell") %>% setkey( round_bp) %>%
    data.table:::split.data.table(by = c("chr", "cell"), flatten = F, keep.by = F) 
  
  smooth_distance_chr_cell.function <- function(chromo){
    
    cell_pairwise <- CJ(as.integer(names(chromo)), 
                        as.integer(names(chromo))
    )[V1 < V2] %>%  split(.$V1)
    
    dists <- map(cell_pairwise, function(x){
      
      data.table:::merge.data.table(chromo[[as.character(x$V1[1])]], 
                                    rbindlist(chromo[as.character(x$V2)], idcol = "cell"), 
                                    c("round_bp")
      )[,.(dens_min = sum(pmin(  density.x, density.y)) / 
             min(c(tot_density.x[1], tot_density.y[1]))), 
        .(cell)]              
    })
    
    return(dists)}
  
  smooth_distance_chr_cell <- map(smooth_chr_cell, smooth_distance_chr_cell.function)
  
  dist_matrix <- rbindlist(map(smooth_distance_chr_cell, rbindlist, idcol = "from"), 
                           idcol = "chr")[chr != "X", .(dist = sum(dens_min) / 22), .(from, cell)] %>% 
    mutate(from = as.integer(from),
           cell = as.integer(cell)) %>% 
    full_join(CJ(from = unique(c(.$from, .$cell)), cell = unique(c(.$from, .$cell))),
              by = c("from", "cell")) %>% 
    mutate(dist = ifelse(from == cell, 1, round(dist, 3)),
           dist = dist^-1 - 1,
           dist = ifelse(dist == Inf | is.na(dist), 1000, dist)) %>%
    arrange(from, cell)
  
  return(dist_matrix)
  
}

super_sphase_order       <- function(inputfile){
  
  dist_matrix <- inputfile %>% 
    pivot_wider(names_from = from, values_from = dist) %>%
    column_to_rownames("cell")
  
  umap.custom                <- umap.defaults
  umap.custom$n_neighbors    =  round(log2(nrow(dist_matrix)))
  umap.custom$n_components   =  2
  umap.custom$n_epochs       =  200
  umap.custom$input          =  "dist"
  umap.custom$verbose        =  F
  
  umap.custom2               <- umap.custom
  umap.custom2$n_neighbors   =  round(log2(nrow(dist_matrix)) * 2)
  umap.custom2$n_components  =  1
  
  umap_matrix <- mclapply(1:100, function(x){
    
    set.seed(1000 + x)
    cell_umap2      <- umap(as.matrix(as.dist(dist_matrix)), config = umap.custom2)
    
    scale(cell_umap2$layout)
    
  }, mc.cores = 5)
  
  bootstrap_threshold <- 0.85
  cluster_threshold   <- 0.8
  
  run_select <- colMeans(cor(bind_cols(umap_matrix), method = "spearman")) > bootstrap_threshold
  
  if(sum(run_select) == 0) { print("better check your bootstrap results")
    return( colMeans(cor(bind_cols(umap_matrix), method = "spearman")))}
  
  umap_tour_bootstrap <- as.data.table(
    data.frame(bind_cols(umap_matrix[run_select]), 
               cell = names(dist_matrix)) %>% 
      pivot_longer(-cell)
  )[order(value), cluster := .(cumsum(c(0, diff(value) > 0.1))), .(cell)
  ][, N := .(.N / sum(run_select)), .(cell, cluster)
  ][, N_threshold := .(max(N) > cluster_threshold), .(cell)
  ][, mean := .(weighted.mean(value, N > ifelse(N_threshold[1], cluster_threshold, 0))), .(cell)
  ][, stdev := .(sd(value[N > cluster_threshold])), .(cell)
  ][, rank := .(as.double(dense_rank(mean))), .(N_threshold)
  ][, rank := .(rank / max(rank)), .(N_threshold)] 
  
  selected_cells <- umap_tour_bootstrap[(N_threshold) & !is.na(stdev)] %>% 
    distinct(cell, .keep_all = T)
  
  cell_umap  <- umap(as.matrix(as.dist(dist_matrix[selected_cells$cell, selected_cells$cell])), 
                     config = umap.custom)
  
  umap_tour_final <- data.frame(cell_umap$knn$indexes,
                                cell_umap$layout) %>% 
    rownames_to_column(   "cell") %>% 
    left_join(selected_cells[,.(cell, mean, stdev, rank, N)], by = "cell") %>%
    pivot_longer(c(-cell, -X1.1, -X2.1, -mean, -stdev, -rank, -N)) %>%
    mutate(X2 = cell_umap$layout[value,1],
           Y2 = cell_umap$layout[value,2])
  
  return(list(bootstrap = umap_tour_bootstrap, 
              final_order = umap_tour_final))
}

super_pc                 <- function(inputfile, maxx = 1e6){
  
  inputfile %>%
    data.table:::split.data.table(by = "cell") %>%
    map( function(single_cell){ single_cell[, dist_1d(bp, maxx), .(chr)]}) %>% 
    return()
  
}

super_fit                <- function(inputfile){
  
  
  better_fit.function <- function(cell_index, minexp = 1000, maxd = 4e5){

    A <- inputfile[[cell_index]][dist < maxd & dist > 0
    ][, .(N2 = as.double(.N + runif(1, 0, 0.01))), keyby = .(dist = round((dist + 2499.9) / 5000) * 5000)
    ][, N := frollmean(N2, 7, align = "center")
    ][order(dist), diff1 := c(0,diff(sign(diff(N))), 0)
    ][,Nscale := (N - min(N, na.rm = T)) / (max(N2, na.rm = T) - min(N, na.rm = T))]
    
    AA <- A[(diff1 == -2)  & Nscale > 0.05
    ][order(-N)
    ][,rep(dist, 20)]
    
    init_fail <- length(AA) == 0
    
    if(length(AA) == 0){ 
      
      AA <- rep(seq(10e4, 25e4, length.out = 20), 20) }
    
    if(nrow(inputfile[[cell_index]]) < 1000){return(NULL)}
    
    rbind(
      map_dfr(1:20, function(x){
        EM_doublemix(inputfile[[cell_index]][dist < maxd & dist > 0
        ][sample(.N, round(.N / 100 * 2^(-log10(.N) + log2(320))))]$dist, 
        c(2e4, maxd, 0.25, 0.7), 
        tol = 1e-8, minexppar = 1000) %>%
          {as.data.table(t(c(.$pars_priors, max(.$loglike[.$loglike != 0]), NA)))}}
      )[, .(parameter = c(1,6,7,10,11),
            Estimate = map_dbl(.SD, mean),
            std_err = map_dbl(.SD, sd),
            stde_mean = NA,
            est_sd = NA, 
            model = "double"), .(peak_init = V6)],
      
      map_dfr(c(AA), function(x){
        EM_triplemix(inputfile[[cell_index]][dist < maxd & dist > 0
        ][sample(.N, round(.N / 100 * 2^(-log10(.N) + log2(320))))]$dist, 
        c(2e4, x + rnorm(1, sd = 1000), (2e4 + init_fail * 5e3)^2, maxd, 0.2, 0.1, 0.7), 
        tol = 1e-8, minexppar = 1000) %>%
          {as.data.table(t(c(.$pars_priors, max(.$loglike[.$loglike != 0]), x)))}}
      )[, .(parameter = c(1,4:7,9:10,11),
            Estimate = map_dbl(.SD, mean),
            std_err = map_dbl(.SD, sd),
            stde_mean = sd(V2),
            est_sd = mean(sqrt(V3)),
            model = "triple"), .(peak_init = V9)],
      
      map_dfr(AA, function(x){
        EM_quadruplemix(inputfile[[cell_index]][dist < maxd & dist > 0
        ][sample(.N, round(.N / 100 * 2^(-log10(.N) + log2(320))))]$dist, 
        c(2e4, 2.5e4 + rnorm(1, sd = 1000), 2e4^2, 
          x + rnorm(1, sd = 1000), (2e4 + init_fail * 5e3)^2, maxd, 0.2, 0.05, 0.15, 0.6), 
        tol = 1e-8, minexppar = 1000) %>%
          {as.data.table(t(c(.$pars_priors, max(.$loglike[.$loglike != 0]), x)))}}
        
      )[, .(parameter = 1:11,
            Estimate = map_dbl(.SD, mean),
            std_err = map_dbl(.SD, sd),
            stde_mean = sd(V4),
            est_sd = mean(sqrt(V5)),
            model = "quadruple"), .(peak_init = V12)])
    
  }
  
  setNames(names(inputfile), names(inputfile)) %>%
    mclapply(better_fit.function, mc.cores = 5) %>%
    rbindlist(idcol = "cell") %>% return()
  
}

super_fit2               <- function(inputfile){
  
  
  better_fit.function <- function(cell_index, minexp = 1000, maxd = 1e6){
    
    A <- inputfile[[cell_index]][dist < maxd & dist > 0
    ][, .(N2 = as.double(.N + runif(1, 0, 0.01))), keyby = .(dist = round((dist + 2499.9) / 5000) * 5000)
    ][, N := frollmean(N2, 7, align = "center")
    ][order(dist), diff1 := c(0,diff(sign(diff(N))), 0)
    ][,Nscale := (N - min(N, na.rm = T)) / (max(N2, na.rm = T) - min(N, na.rm = T))]
    
    AA <- A[(diff1 == -2)  & Nscale > 0.05
    ][order(-N)
    ][,rep(dist, 20)]
    
    init_fail <- length(AA) == 0
    
    if(length(AA) == 0){ 
      
      AA <- rep(seq(10e4, 25e4, length.out = 20), 20) }
    
    if(nrow(inputfile[[cell_index]]) < 1000){return(NULL)}
    
    rbind(
      map_dfr(1:20, function(x){
        EM_doublemix(inputfile[[cell_index]][dist < maxd & dist > 0
        ][sample(.N, round(.N / 100 * 2^(-log10(.N) + log2(320))))]$dist, 
        c(2e4, maxd, 0.25, 0.7), 
        tol = 1e-8, minexppar = 1000) %>%
          {as.data.table(t(c(.$pars_priors, max(.$loglike[.$loglike != 0]), NA)))}}
      )[, .(parameter = c(1,6,7,10,11),
            Estimate = map_dbl(.SD, mean),
            std_err = map_dbl(.SD, sd),
            stde_mean = NA,
            est_sd = NA, 
            model = "double"), .(peak_init = V6)],
      
      map_dfr(c(AA), function(x){
        EM_triplemixf(inputfile[[cell_index]][dist < maxd & dist > 0
        ][sample(.N, round(.N / 100 * 2^(-log10(.N) + log2(320))))]$dist, 
        c(2e4, 0, 3e5^2, maxd, 0.2, 0.1, 0.7), 
        tol = 1e-8, minexppar = 1000) %>%
          {as.data.table(t(c(.$pars_priors, max(.$loglike[.$loglike != 0]), x)))}}
      )[, .(parameter = c(1,2,3,6:8,10,11),
            Estimate = map_dbl(.SD, mean),
            std_err = map_dbl(.SD, sd),
            stde_mean = sd(V2),
            est_sd = mean(sqrt(V3)),
            model = "triple"), .(peak_init = V9)],
      
      map_dfr(AA, function(x){
        EM_quadruplemixf(inputfile[[cell_index]][dist < maxd & dist > 0
        ][sample(.N, round(.N / 100 * 2^(-log10(.N) + log2(320))))]$dist, 
        c(2e4, 0, 3e5^2, x + rnorm(1, sd = 1000), (2e4 + init_fail * 5e3)^2, maxd, 0.15, 0.3, 0.15, 0.4), 
        tol = 1e-8, minexppar = 1000) %>%
          {as.data.table(t(c(.$pars_priors, max(.$loglike[.$loglike != 0]), x)))}}
        
      )[, .(parameter = 1:11,
            Estimate = map_dbl(.SD, mean),
            std_err = map_dbl(.SD, sd),
            stde_mean = sd(V4),
            est_sd = mean(sqrt(V5)),
            model = "quadruple"), .(peak_init = V12)])
    
  }
  
  setNames(names(inputfile), names(inputfile)) %>%
    mclapply(better_fit.function, mc.cores = 5) %>%
    rbindlist(idcol = "cell") %>% return()
  
}

super_chr_bs_fit         <- function(inputfile){
  
  inputfits  <- fread(paste0("../pc_correct/",   sub("_.*", "_pc_correct.tsv",   sub("../[[:graph:]_]{1,20}/", "", arguments[2]))))[(diffusion)]
  inputfits2 <- fread(paste0("../model_select/", sub("_.*", "_model_select.tsv.gz", sub("../[[:graph:]_]{1,20}/", "", arguments[2]))))
  
  better_fit.function <- function(cell_index, minexp = 1000, maxd = 1e6){
  
    init_peak <- inputfits[cell == cell_index]
    init_peak2 <- inputfits2[init_peak, on = c("diffusion", "cell", "peak_init", "model")][parameter != 11]
    
    if(nrow(inputfile[[cell_index]]) < 1000){return(NULL)}

    listout <- try({
      list(
        (inputfile[[cell_index]] %>%
           split(.$chr) %>%
           map_dfr(function(x){
             EM_doublemix(x[dist < maxd & dist > 0]$dist, 
                          c(2e4, maxd, 0.25, 0.7), tol = 1e-7, minexppar = 1000) %>%
               {as.data.table(t(c(.$pars_priors, max(.$loglike[.$loglike != 0])))
               )[,chr := x$chr[1]
               ][, n_obs := nrow(x)][]}}
           ))[,.(chr, V1, V6 = V2, V7 = V3, V10 = V4,
                 loglike = V5, model = "double", peak_init = init_peak$Estimate, n_obs)],
        
        (inputfile[[cell_index]] %>%
           split(.$chr) %>%
           map_dfr(function(x){
             EM_triplemixf(x[dist < maxd & dist > 0]$dist, 
                           c(2e4, 0, 3e5^2, maxd, 0.2, 0.1, 0.7), 
                           tol = 1e-7, minexppar = 1000) %>%
               {as.data.table(t(c(.$pars_priors, max(.$loglike[.$loglike != 0]))))[,chr := x$chr[1]][, n_obs := nrow(x)][]}}
           ))[,.(chr, V1, V2, V3, V6 = V4, V7 = V5, V8 = V6, V10 = V7,
                 loglike = V8, model = "triple", peak_init = init_peak$Estimate, n_obs)],
        
        (inputfile[[cell_index]] %>%
           split(.$chr) %>%
           map_dfr(function(x){
             EM_quadruplemixf(x[dist < maxd & dist > 0]$dist, init_peak2$Estimate, 
                              tol = 1e-7, minexppar = 1000) %>%
               {as.data.table(t(c(.$pars_priors, max(.$loglike[.$loglike != 0]))))[,chr := x$chr[1]][, n_obs := nrow(x)][]}}
           ))[,.(chr, V1, V2, V3, V4, V5, V6, V7, V8, V9, V10,
                 loglike = V11, model = "quadruple", peak_init = init_peak$Estimate, n_obs)]
        
      ) %>% rbindlist(fill = T)
    })
    
    gc()
    
    return(listout)
    
  }
  
  chr_fits <- setNames(names(inputfile), names(inputfile)) %>%
    
    {.[. %in% inputfits$cell]} %>%
    
    mclapply(better_fit.function, mc.cores = 3) %>% 
    rbindlist(idcol = "cell") 
  
  chr_fits[, BIC := setNames(c(2, 5, 8), c("double", "triple", "quadruple"))[model] * log(n_obs) - loglike[1], .(cell, chr, model)]
  chr_fits[, BICrank := frank(BIC, ties.method = "dense"), .(cell, chr)]
  
  return(chr_fits)
  
}

super_model_select       <- function(inputfile){
  
  rbind(
    fread(paste0("../better_fit/", sub("_.*", "_better_fit.tsv", sub("../[[:graph:]_]{1,20}/", "", arguments[2])))
    )[parameter != 11, loglike := EM_loglike(inputfile[[as.character(cell)]][dist > 0 & dist < 4e5]$dist, Estimate), 
      .(cell, model, peak_init)
    ][, n_obs := nrow(inputfile[[as.character(cell)]]), .(cell, model, peak_init)
    ][, BIC := setNames(c(2, 5, 8), c("double", "triple", "quadruple"))[model] * log(n_obs) - loglike[1], 
      .(cell, model, peak_init)
    ][,diffusion := F],
    fread(paste0("../better_fit2/", sub("_.*", "_better_fit2.tsv", sub("../[[:graph:]_]{1,20}/", "", arguments[2])))
    )[parameter != 11, loglike := EM_loglikef(inputfile[[as.character(cell)]][dist > 0 & dist < 1e6]$dist, Estimate), 
      .(cell, model, peak_init)
    ][, n_obs := nrow(inputfile[[as.character(cell)]]), .(cell, model, peak_init)
    ][, BIC := setNames(c(2, 4, 7), c("double", "triple", "quadruple"))[model] * log(n_obs) - loglike[1], 
      .(cell, model, peak_init)
    ][,diffusion := T]
  )
  
}

super_hmm_fit            <- function(inputfile){
  
  exp_hmm_cell.function <- function(single_cell){
    
    sc <- single_cell[order(chr, bp)][, dist := as.integer(c(0, diff(bp))), .(chr)]
    
    sc_response <- flexmix(dist ~ 1, model = list(FLXMCdist1(dist = "exp")), k = 2, 
                           data = sc[dist > 0])
    
    if(length(parameters(sc_response)) < 2){return(NULL)}
    
    
    sc_hsmm_model <- hsmmspec(init = c(0.5, 0.5),
                              transition = matrix(c(0, 1, 1, 0), 2), 
                              parms.emission = list(lambda = 1 / sort(parameters(sc_response))),
                              sojourn  = list(d = cbind(dunif(1:200, 1, 50),
                                                        dunif(1:200, 1, 200)),
                                              type = "gamma"),
                              dens.emission = dexp.hsmm, 
                              rand.emission = rexp.hsmm, 
                              mstep = mstep.exp) 
    
    sc_hsmm_data <- list(x = sc[dist > 0]$dist,
                         N = sc[dist > 0][,.(N = .N), .(chr)]$N)
    class(sc_hsmm_data) <- "hsmm.data"
    
    sc_hsmm_fit <- NULL
    
    try(
      sc_hsmm_fit <- hsmmfit(x     = sc_hsmm_data, 
                             model = sc_hsmm_model, 
                             mstep = mstep.exp,
                             graphical = F,
                             M = 200))
    
    if(class(sc_hsmm_fit) != "hsmm"){
      
      sc_hsmm_model <- hsmmspec(init = c(0.5, 0.5),
                                transition = matrix(c(0, 1, 1, 0), 2), 
                                parms.emission = list(lambda = 1 / sort(parameters(sc_response))), 
                                sojourn  = list(d = cbind(dunif(1:300, 1, 50),
                                                          dunif(1:300, 1, 300)),
                                                type = "gamma"),
                                dens.emission = dexp.hsmm, 
                                rand.emission = rexp.hsmm, 
                                mstep = mstep.exp) 
      
      
      try(
        sc_hsmm_fit <- hsmmfit(x     = sc_hsmm_data, 
                               model = sc_hsmm_model, 
                               mstep = mstep.exp,
                               graphical = F,
                               M = 300))}
    
    if(class(sc_hsmm_fit) != "hsmm"){return(class(sc_hsmm_fit))}
    
    sc_hsmm_predict <- predict(sc_hsmm_fit, sc_hsmm_data, "smoothed")  
    
    list(data = data.table:::merge.data.table(sc, by = c("chr", "bp"), all.x = T,
                                              cbind(sc[dist > 0][,.(chr, bp)], 
                                                    data.table(state     = sc_hsmm_predict$s, 
                                                               posterior = sc_hsmm_predict$p[,2]))),
         model = unlist(list(init = sc_hsmm_fit$model$init, 
                             emission = sc_hsmm_fit$model$parms.emission, 
                             sojourn_shape = sc_hsmm_fit$model$sojourn$shape,
                             sojourn_scale = sc_hsmm_fit$model$sojourn$scale)))
    
  }
  
  inputfile[,.(cell, chr, bp, allele)] %>%
    data.table:::split.data.table(by = "cell", keep.by = F) %>%
    mclapply(exp_hmm_cell.function, mc.cores = 10) %>% return() #mc.cores = 10
  
}

super_single_fork        <- function(inputfile){
  
  inputfile[order(cell, chr, bp, state)
  ][, state_after := lead(state), .(cell, chr)
  ][is.na(state), state := state_after
  ][is.na(state_after), state_after := state  
  ][, group := cumsum((state == 1 & state_after == 2) | 
                        (state == 2 & state_after == 1)), .(cell, chr)
  ][state == 2 & state_after == 1, group := group - 1
  ][,.(#bp    = list(as.integer(bp)), 
    #dist  = list(as.integer(dist)[-1]), 
    state = max(state),
    allele_asigned = case_when(!any(!allele == "none") ~ "none",
                               ("A" %in% allele & "B" %in% allele) | ("both" %in% allele) ~ "both",
                               "A" %in% allele ~ "A",
                               "B" %in% allele ~ "B"),
    N = .N,
    center = min(bp) + sum(dist[-1]) / 2,
    width = sum(dist[-1])),
    .(cell, chr, group)]
  
}

# simulated pair correlation

pc_sim <- CJ(U = seq(0.5, 3, by = 0.05), 
            S = seq(0.1, 0.9, by = 0.02)
  )[,.(U = rep(U, 10), S = rep(S, 10))
  ][, i := 1:.N
  ][, c("u", "s") := pc_sim(c(U, S), 2000), .(i)]

pc_sim[, c("u_R", "s_R")       := .(abs(u - mean(u)) / sd(u), abs(s - mean(s)) / sd(s)), .(U,S)]
pc_sim[, c("u_rank", "s_rank") := .(rank(u_R), rank(s_R)), .(U,S)]

U_model <- loess(U ~ u + s, pc_sim[u_rank < 10 | s_rank < 10])
S_model <- loess(S ~ u + s, pc_sim[u_rank < 10 | s_rank < 10])

all_fits <- list.files("hpc/vincentvb/forkdot/model_select/", ".tsv.gz") %>%
  map(function(x){
    print(x)
    
    fread(paste0("hpc/vincentvb/forkdot/model_select/", x)
    )[((!diffusion) & model != "quadruple") | (diffusion)
    ][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
    ][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(cell, diffusion, model, peak_init)
    ][(failed), .SD[BIC == min(BIC)], .(diffusion, cell)
    ][parameter == 4][, exp := x]
  }) %>% rbindlist()

all_fits[, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))]
all_fits[, s := predict(S_model, data.table(u = Estimate / 75000, s = est_sd / 75000))]

# expressed gene coverage marginal
hg38.gtf    <- fread("hpc/vincentvb/forkdot/Homo_sapiens.GRCh38.106.chr.gtf.gz")

hg38_gene_location <- setNames(hg38.gtf[,c(1,4,5,9)], c("chr", "from", "to", "info")
)[,.(chr, from, to, geneid = str_extract(info, "gene_id \"ENSG[0-9]{11}\""),  
     biotype = str_extract(info, "\"protein_coding\""))
][, .(GENEID = gsub("gene_id |\"", "", geneid[1]), 
      range = max(to - from)), .(chr, from, to, biotype)
][, .(from = from[which.max(range)], to = to[which.max(range)]), .(GENEID, chr, biotype)
][biotype == "\"protein_coding\""
][order(chr, from)
][,overlap := lag(to) - from, .(chr) #to - lead(from), na.rm = T)
][is.na(overlap), overlap := 0
][, group := cumsum(overlap < 100)
][,n := .N,.(chr, group)]

scEU <- fread("hpc/vincentvb/forkdot/scEU_pulse60min_combined_V2.tsv.gz"
)[, cell_cycle := fread("hpc/vincentvb/forkdot/scEU_index.tsv.gz"
)[,setNames(paste0(Cell_cycle_relativePos > 0.333,  Cell_cycle_relativePos < 0.75), Cell_Id)
][cell]
][,.(count = sum(count)), .(cell_cycle, Gene_Id, sample)
][cell_cycle == "TRUETRUE" & sample == "nascent",.(GENEID = Gene_Id, value = count)]

gene_loc <- merge.data.table(hg38_gene_location, scEU, by = "GENEID", sort = F
)[order(chr, from)
][,overlap := lag(to) - from, .(chr) #to - lead(from), na.rm = T)
][is.na(overlap), overlap := 0
][, group := cumsum(overlap < 100)
][,c("n", "order") := .(.N, 1:.N), .(chr, group)
][, c("from2", "to2") := .(fifelse(value < lag(value)  & order > 1,  lag(to), from),
                           fifelse(value < lead(value) & order < .N, lead(from), to)), .(chr, group)
][(to2 - from2) > 10
][,.(chr, GENEID, from = from2, to = to2, 
     value = round(log10(value + 1))) #, d = 0
] %>% setkeyv(c("chr", "from", "to"))

single_forks <- fread("hpc/vincentvb/forkdot/single_fork/JvB-exp049-15m-EdU-beads_single_fork.tsv.gz")

sphase <- as.data.table(readRDS("hpc/vincentvb/forkdot/sphase_order/JvB-exp049-15m-EdU-beads_sphase.RDS"
)$final_order)[, unique(rank), .(cell)][,setNames(V1, cell)]

coverage_p_cell <- foverlaps(
  single_forks[state == 2
  ][, width := width + width / N
  ][,.(chr, from = center - width / 2, to = center + width / 2, width, 
       tot_covered = sum(width)), .(cell)],
  gene_loc, 
  by.x = c("chr", "from", "to"), 
  by.y = c("chr", "from", "to")
)[!is.na(value),.(gene_overlap = sum(pmin(i.to, to) - pmax(i.from, from))), .(cell, tot_covered, value) 
][, rank := sphase[as.character(cell)]]


{
  A <- scEU[, .(N = .N), .(group = round(log10(value + 1)), value = round(log10(value + 1) * 10) / 10)
  ][group < 4] %>%
    ggplot(aes(x = value, fill = as.factor(group), y = N)) + geom_col() + 
    theme_bw() +
    theme(legend.position = "none", plot.background = element_blank())
  
  B <- coverage_p_cell[value < 4
  ][order(rank), gene_overlap := frollmean(gene_overlap / sum(gene_overlap), 25, align = "center"),.(value)
  ]%>%
    ggplot() + 
    geom_line(aes(x = rank, y = gene_overlap, col = as.factor(value))) +
    theme_bw() +
    theme(legend.position = "none") 
  
  C <- coverage_p_cell[value < 4] %>%
    ggplot() + 
    geom_col(aes(x = rank, y = gene_overlap / tot_covered, fill = as.factor(value))) +
    theme_bw() +
    theme(legend.position = "bottom") 
  
  
  png();B;dev.off()
  
}












