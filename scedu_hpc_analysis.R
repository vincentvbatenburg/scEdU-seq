#!/usr/bin/env Rscript


cat("\n#########################################
     \nSTART START START START START START START
     \n#########################################\n")


# packages
suppressPackageStartupMessages(library(data.table))         
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(slurmR))
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(flexmix))
suppressPackageStartupMessages(library(mhsmm))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(umap))

# set package options
setDTthreads(threads = 0, restore_after_fork = T)


# dependencies

path2source = "/hpc/hub_oudenaarden/vincentvb/forkdotV2/src/"

#sourceCpp("/hpc/hub_oudenaarden/vincentvb/forkdotV2/src/EM_triplemix.cpp")
#sourceCpp("/hpc/hub_oudenaarden/vincentvb/forkdotV2/src/dist_1d.cpp")
#sourceCpp("/hpc/hub_oudenaarden/vincentvb/forkdotV2/src/Utility.cpp")
#sourceCpp("/hpc/hub_oudenaarden/vincentvb/forkdotV2/src/pc_correct.cpp")

Slurm_EvalQ_custom <- function(expr, njobs = 2L, job_name = opts_slurmR$get_job_name(), 
                               tmp_path = opts_slurmR$get_tmp_path(), plan = "collect", 
                               sbatch_opt = list(), rscript_opt = list(), seeds = NULL, 
                               compress = TRUE, export = NULL, export_env = NULL, libPaths = .libPaths(), 
                               hooks = NULL, overwrite = TRUE, preamble = NULL) 
{
  plan <- the_plan(plan)
  # added the try to silence the warning of temp path not existing
  try({check_full_path(tmp_path = tmp_path, job_name = job_name, 
                  overwrite = overwrite)})
  sbatch_opt <- check_sbatch_opt(sbatch_opt, job_name = job_name, 
                                 ntasks = 1L)
  opts_slurmR$set_tmp_path(tmp_path)
  opts_slurmR$set_job_name(job_name)
  sexpr <- deparse(substitute(expr))
  if (is.null(export_env)) 
    export_env <- parent.frame()
  rscript <- new_rscript(njobs, libPaths = libPaths, tmp_path = tmp_path, job_name = job_name, pkgs = c(list_loaded_pkgs(), parallel = "/hpc/hub_oudenaarden/vincentvb/bin/miniconda3/envs/forkdotV2/lib/R/library"))
  if (length(export)) {
    rscript$add_rds(mget(export, envir = export_env), compress = compress) # removed the index = TRUE as the add_rds function doesnt have an index argument
  }
  rscript$set_seed(seeds)
  rscript$append(paste0("ans <- {list(\n", paste0(gsub("^", 
                                                       "   ", sexpr), collapse = "\n"), "\n)}"))
  rscript$finalize("ans", compress = compress)
  rscript$write()
  bash <- new_bash(njobs = njobs, job_name = job_name, output = snames("out", 
                                                                       job_name = job_name, tmp_path = tmp_path), filename = snames("sh", 
                                                                                                                                    job_name = job_name, tmp_path = tmp_path))
  bash$add_SBATCH(sbatch_opt)
  bash$append(c(opts_slurmR$get_preamble(), preamble))
  bash$Rscript(file = snames("r", job_name = job_name, tmp_path = tmp_path), 
               flags = rscript_opt)
  bash$write()
  hooks <- c(hooks, list(function(res, job, ...) {
    stats::setNames(res, paste0(job$jobid, "_", 1:job$njobs))
  }))
  ans <- new_slurm_job(call = match.call(), rscript = snames("r", 
                                                             job_name = job_name, tmp_path = tmp_path), bashfile = snames("sh", 
                                                                                                                          job_name = job_name, tmp_path = tmp_path), robjects = NULL, 
                       njobs = njobs, opts_job = sbatch_opt, opts_r = opts_slurmR$get_opts_r(), 
                       hooks = hooks)
  if (plan$collect) 
    return(Slurm_collect(sbatch(ans, wait = plan$wait, submit = plan$submit)))
  else return(sbatch(ans, wait = plan$wait, submit = plan$submit))
}

environment(Slurm_EvalQ_custom) <- asNamespace('slurmR')

mstep.exp <- function(x, wt)       {list(lambda = apply(wt, 2, function(w){weighted.mean(x, w)}))}
rexp.hsmm <- function(j, model)    {rexp(1, 1 / model$parms.emission$lambda[j])}
dexp.hsmm <- function(x, j, model) {dexp(x, 1 / model$parms.emission$lambda[j])}

# dealing with script arguments
arguments <- commandArgs(trailingOnly=TRUE)

print(arguments)

jobid  = gsub(".tsv.gz|pl[0-9]{1,2}", "", arguments[1])

config <- NULL

try({config <- fread(arguments[2], header = F) }, silent = T)

print(config)

# making new folder for all pipelines to run in
if(!str_detect(getwd(), jobid)){
  
  if(!any(grepl(jobid, list.dirs()))) system(paste("mkdir", jobid))
  
  for(file in list.files(pattern = jobid) %>% {.[str_detect(., "tsv.gz")]}) system(paste("mv", file, jobid))
  
  setwd(jobid)
  
  opts_slurmR$set_tmp_path(getwd())
  
}

# checking which files already exist and which parts to run
file.checker <- function(){
  c("poissontest", "filterer", "overlapscore", "sphaseorder", "pc", "psmashc", "mmfit", "BICfit", "hmmfit", "singlefork") %>%
    {setNames(str_detect(paste0(list.files(), collapse = "|"), paste0("_", .)), .)}
}

file_check <- file.checker()

print(file_check)

# pipeline functions
super_poissontest  <- function(exp, species = "human", path2blacklist = NULL, is_allelic = TRUE){
  
  setDTthreads(threads = 0, restore_after_fork = T)
  
  try({for(i in 1:nrow(config)) assign(config$V1[i], config$V2[i])}, silent = T)
  
  totalbins <- c(human = 29111, mouse = 25627)[species]
  
  if(!is.null(path2blacklist)) blacklistbins <- fread(path2blacklist)
  
  cuts_cell <- rbindlist(map(list.files(pattern = exp) %>% {.[str_detect(., ".tsv.gz") & str_detect(., "-pl[0-9]{1,2}")]}, 
                             function(x){
                               plate <- fread(x)
                               
                               if(is_allelic){
                                 
                                 if(ncol(plate) != 6){ 
                                   write_tsv(data.table(error = paste0(x, "__not_6_columns")), file = "error.tsv")
                                   stop()} 
                                 
                                 plate_pros <- setNames(plate, c("cell", "chr", "bp", "strand", "readlength", "allele")
                                 )[chr %in% c(1:22, "X") & readlength != "None" & as.numeric(readlength) < 1e3]
                               } else {
                                 
                                 if(ncol(plate) != 5){ 
                                   write_tsv(data.table(error = paste0(x, "__not_5_columns")), file = "error.tsv")
                                   stop()} 
                                 
                                 plate_pros <- setNames(plate, c("cell", "chr", "bp", "strand", "readlength")
                                 )[chr %in% c(1:22, "X") & readlength != "None" & as.numeric(readlength) < 1e3]
                               }
                               
                               if(is.null(path2blacklist)){
                                 return(plate_pros)
                               } else { 
                                 return(plate_pros[,.SD[!((round(bp / 5000) * 5000)  %in% blacklistbins[chrbb == chr]$bp_round)], .(chr)])
                               }
                             })
  )[, bp_round := as.integer(round((bp + 49999.99) / 100000) * 100000 - 49999.99)
  ][,.(N = as.double(.N)), .(cell, chr, bp_round)
  ][,.(mean = sum(N) / totalbins,
       variance = sum((N - (sum(N) / totalbins))^2) / totalbins), .(cell)
  ][, poisson_score := log(sqrt(variance) / mean) + 0.5 * log(mean)]
  
  write_tsv(cuts_cell, paste0(exp, "_poissontest.tsv.gz"))
  
}

super_filterer     <- function(exp, species = "human", cpus = 4, poisson_threshold = 0.1, mean_lower = -1, mean_upper = 2.5, path2blacklist = NULL, is_allelic = TRUE){
  
  setDTthreads(threads = 0, restore_after_fork = T)
  
  try({for(i in 1:nrow(config)) assign(config$V1[i], config$V2[i])}, silent = T)
  
  if(!is.null(path2blacklist)) blacklistbins <- fread(path2blacklist)
  
  cuts_cell <- rbindlist(map(list.files(pattern = exp) %>% {.[str_detect(., ".tsv.gz") & str_detect(., "-pl[0-9]{1,2}")]}, 
                             function(x){
                               plate <- fread(x)
                               
                               if(is_allelic){
                                 
                                 if(ncol(plate) != 6){ 
                                   write_tsv(data.table(error = paste0(x, "__not_6_columns")), file = "error.tsv")
                                   stop()} 
                                 
                                 plate_pros <- setNames(plate, c("cell", "chr", "bp", "strand", "readlength", "allele")
                                 )[chr %in% c(1:22, "X") & readlength != "None" & as.numeric(readlength) < 1e3]
                               } else {
                                 
                                 if(ncol(plate) != 5){ 
                                   write_tsv(data.table(error = paste0(x, "__not_5_columns")), file = "error.tsv")
                                   stop()} 
                                 
                                 plate_pros <- setNames(plate, c("cell", "chr", "bp", "strand", "readlength")
                                 )[chr %in% c(1:22, "X") & readlength != "None" & as.numeric(readlength) < 1e3]
                               }
                               
                               if(is.null(path2blacklist)){
                                 return(plate_pros)
                               } else { 
                                 return(plate_pros[,.SD[!((round(bp / 5000) * 5000)  %in% blacklistbins[chrbb == chr]$bp_round)], .(chr)])
                               }
                             })
  )
  
  cuts_cell[,readlength := as.numeric(readlength)]
  
  poisson_test <- fread(list.files(pattern = "_poissontest")
  )[order(cell),cell2 := 1:.N
  ][poisson_score > as.numeric(poisson_threshold) & log(mean) > as.numeric(mean_lower) & log(mean) < as.numeric(mean_upper)]
  
  if(is_allelic){
    exp_fit_cell.function <- function(single_cell){
      
      sc  <- single_cell[,dist := as.integer(c(0, diff(bp))), .(chr)]
      
      if(nrow(sc[dist > 0]) < 10) return(NULL)
      
      YYY <- flexmix(dist ~ 1, model = list(FLXMCdist1(dist = "exp")), k = 2, 
                     data = sc[dist > 0])
      
      merge.data.table(sc[,.(chr, bp, allele, readlength, strand)], 
                       cbind(sc[dist > 0][,.(bp, chr)], 
                             YYY@posterior$scaled[,which.max(parameters(YYY))]),
                       all.x = T, by = c("chr", "bp")
      )[order(bp), posterior := pmax(V2, shift(V2, -1), na.rm = T), .(chr)
      ][, V2 := NULL]
      
    }
  } else {
    exp_fit_cell.function <- function(single_cell){
      
      sc  <- single_cell[,dist := as.integer(c(0, diff(bp))), .(chr)]
      
      if(nrow(sc[dist > 0]) < 10) return(NULL)
      
      YYY <- flexmix(dist ~ 1, model = list(FLXMCdist1(dist = "exp")), k = 2, 
                     data = sc[dist > 0])
      
      merge.data.table(sc[,.(chr, bp, readlength, strand)], 
                       cbind(sc[dist > 0][,.(bp, chr)], 
                             YYY@posterior$scaled[,which.max(parameters(YYY))]),
                       all.x = T, by = c("chr", "bp")
      )[order(bp), posterior := pmax(V2, shift(V2, -1), na.rm = T), .(chr)
      ][, V2 := NULL]
      
    }
  }
  
  exp_fit_cell <- cuts_cell[(cell %in% poisson_test$cell) & (readlength < 1000)
  ][, cell := setNames(poisson_test$cell2, poisson_test$cell)[cell]]  %>% 
    split(.$cell) %>%
    mclapply(exp_fit_cell.function, mc.cores = cpus) %>%
    # {.[map_lgl(.,is.list)]} %>% 
    rbindlist(idcol = "cell")
  
  exp_fit_cell[,posterior := round(posterior * 1000)]
  
  write_tsv(exp_fit_cell, paste0(exp, "_filterer.tsv.gz"))
  
}

super_overlapscore <- function(exp, species = "human"){
  
  setDTthreads(threads = 0, restore_after_fork = T)
  
  try({for(i in 1:nrow(config)) assign(config$V1[i], config$V2[i])}, silent = T)
  
  chrcount <- c(human = 22, mouse = 18)[species]
  
  smooth_chr_cell <- fread(list.files(pattern = "_filterer.tsv.gz")  
  )[posterior > 500, .(cell = as.character(cell), chr, bp)] %>% 
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
                        as.integer(names(chromo)))[V1 < V2] %>% split(.$V1)
    
    dists <- map(cell_pairwise, function(x){
      
      data.table:::merge.data.table(chromo[[as.character(x$V1[1])]], 
                                    rbindlist(chromo[as.character(x$V2)], idcol = "cell"), 
                                    c("round_bp")
      )[,.(dens_min = sum(pmin(  density.x, density.y)) / min(c(tot_density.x[1], tot_density.y[1]))), 
        .(cell)]              
    })
    
    return(dists)}
  
  smooth_distance_chr_cell <- map(smooth_chr_cell, smooth_distance_chr_cell.function)
  
  dist_matrix <- rbindlist(map(smooth_distance_chr_cell, rbindlist, idcol = "from"), idcol = "chr"
  )[chr != "X", .(dist = sum(dens_min) / chrcount), .(from = as.integer(from), cell = as.integer(cell))
  ][CJ(from = unique(c(from, cell)), cell = unique(c(from, cell))),
    on = c("from", "cell")
  ][, dist := round(dist, 3)^-1 - 1
  ][from == cell, dist := 0
  ][is.infinite(dist) | is.na(dist), dist := 1000
  ][order(from, cell)]
  
  write_tsv(dist_matrix, file = paste0(exp, "_overlapscore.tsv.gz"))
  
}

super_sphaseorder  <- function(exp, cluster_threshold = 0.8){
  
  setDTthreads(threads = 0, restore_after_fork = T)
  
  try({for(i in 1:nrow(config)) assign(config$V1[i], config$V2[i])}, silent = T)
  
  dist_matrix <- fread(list.files(pattern = "_overlapscore.tsv.gz")) %>% 
    pivot_wider(names_from = from, values_from = dist) %>%
    column_to_rownames("cell")
  
  umap.custom                <- umap.defaults
  umap.custom$n_neighbors    =  round(log2(nrow(dist_matrix)))
  umap.custom$n_components   =  2
  umap.custom$n_epochs       =  200
  umap.custom$input          =  "dist"
  #  umap.custom$min_dist       =  0.1
  umap.custom$verbose        =  F
  
  umap.custom2               <- umap.custom
  umap.custom2$n_neighbors   =  round(log2(nrow(dist_matrix)) * 2)
  umap.custom2$n_components  =  1
  
  umap_matrix <- mclapply(1:100, function(x){
    
    set.seed(1000 + x)
    cell_umap2      <- umap(as.matrix(as.dist(dist_matrix)), config = umap.custom2)
    
    scale(cell_umap2$layout)
    
  }, mc.cores = 5)
  
  run_select <- colMeans(cor(bind_cols(umap_matrix), method = "spearman")) %>% {. > (max(.) * 0.9)}
  
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
  
  saveRDS(list(bootstrap = umap_tour_bootstrap, 
               final_order = umap_tour_final),
          paste0(exp, "_sphaseorder.RDS"))
}


super_pc           <- function(exp, maxx = 1e6){
  
  setDTthreads(threads = 0, restore_after_fork = T)
  
  try({for(i in 1:nrow(config)) assign(config$V1[i], config$V2[i])}, silent = T)
  
  sourceCpp(paste0(path2source, "dist_1d.cpp"))
  
  fread(list.files(pattern = "_filterer.tsv.gz"))  %>%
    data.table:::split.data.table(by = "cell") %>%
    map( function(single_cell) single_cell[, dist_1d(bp, maxx), .(chr)]) %>% 
    saveRDS(file = paste0(exp, "_pc.RDS"))
  
}

super_psmashc      <- function(exp, bin_size = 5000){
  
  setDTthreads(threads = 0, restore_after_fork = T)
  
  try({for(i in 1:nrow(config)) assign(config$V1[i], config$V2[i])}, silent = T)
  
  readRDS(list.files(pattern = "_pc.RDS")) %>% 
    map(function(x){x[,.(N = as.double(.N)), .(dist = round((dist + (bin_size / 2 - 0.1)) / bin_size) * bin_size)]
    }) %>% 
    rbindlist(idcol = "cell") %>%
    write_tsv(file = paste0(exp, "_psmashc.tsv.gz"))
  
}

super_mmfit        <- function(exp, species = "human", cpus = 4){
  
  try({for(i in 1:nrow(config)) assign(config$V1[i], config$V2[i])}, silent = T)
  
  sourceCpp(paste0(path2source, "EM_triplemix.cpp"))
  
  better_fit.function <- function(celll, minexp = 1000, maxd = 1e6){
    #print(cell_index)
    #sc <- pc_cell[[1]][[2]]
    
    #exp_index = 1; cell_index = "1993"
    
    
    #parameter reminder 
    # double 
    # 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    # em NA NA NA NA um pe NA NA pu
    
    # triple
    # 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    # em NA NA nm nv um pe NA pn pu
    
    # quadruple 
    # 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    # em nm nv nm nv um pe pn pn pu
    
    A <- celll[dist < maxd & dist > 0
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
    
    if(nrow(celll) < 1000){return(NULL)}
    #  print(AA)
    
    rbind(
      map_dfr(1:20, function(x){
        EM_doublemix(celll[dist < maxd & dist > 0
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
        EM_triplemixf(celll[dist < maxd & dist > 0
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
        EM_quadruplemixf(celll[dist < maxd & dist > 0
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
  
  pc_file <- readRDS(list.files(pattern = "_pc.RDS")) #%>%
  
  mclapply(pc_file, better_fit.function, mc.cores = cpus) %>%
    # saveRDS(file = paste0(exp, "_mmfit.rds"))
    rbindlist(idcol = "cell") %>% 
    write_tsv(file = paste0(exp, "_mmfit.tsv.gz"))
  
  
} #number 2

super_BICfit       <- function(exp){
  
  setDTthreads(threads = 0, restore_after_fork = T)
  
  try({for(i in 1:nrow(config)) assign(config$V1[i], config$V2[i])}, silent = T)
  
  sourceCpp(paste0(path2source, "EM_triplemix.cpp"))
  
  pc_input <- readRDS(list.files(pattern = "_pc.RDS"))
  
  fread(list.files(pattern = "_mmfit.tsv.gz")
  )[parameter != 11, loglike := EM_loglikef(pc_input[[as.character(cell)]][dist > 0 & dist < 1e6]$dist, Estimate), 
    .(cell, model, peak_init)
  ][, n_obs := nrow(pc_input[[as.character(cell)]]), .(cell, model, peak_init)
  ][, BIC := setNames(c(2, 4, 7), c("double", "triple", "quadruple"))[model] * log(n_obs) - loglike[1], 
    .(cell, model, peak_init)
  ][,diffusion := T] %>%
    write_tsv(file = paste0(exp, "_BICfit.tsv.gz"))
}


super_hmmfit       <- function(exp, cpus = 4, is_allelic = TRUE){
  
  setDTthreads(threads = 0, restore_after_fork = T)
  
  try({for(i in 1:nrow(config)) assign(config$V1[i], config$V2[i])}, silent = T)
  
  mstep.exp <- function(x, wt)       {list(lambda = apply(wt, 2, function(w){weighted.mean(x, w)}))}
  rexp.hsmm <- function(j, model)    {rexp(1, 1 / model$parms.emission$lambda[j])}
  dexp.hsmm <- function(x, j, model) {dexp(x, 1 / model$parms.emission$lambda[j])} 
  
  if(is_allelic){
  exp_hmm_cell.function <- function(single_cell){
    
    sc <- single_cell[,.(chr, start = fifelse(strand == 1, bp, bp - readlength), 
                         end   = fifelse(strand == 1, bp + readlength, bp), allele) 
    ][order(chr, start * 0.5 + end * 0.5)
    ][, dist := c(0, diff(start * 0.5 + end * 0.5)), .(chr)]
    
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
    
    list(data = data.table:::merge.data.table(sc, by = c("chr", "start", "end"), all.x = T,
                                              cbind(sc[dist > 0][,.(chr, start, end)], 
                                                    data.table(state     = sc_hsmm_predict$s, 
                                                               posterior = sc_hsmm_predict$p[,2]))),
         model = unlist(list(init = sc_hsmm_fit$model$init, 
                             emission = sc_hsmm_fit$model$parms.emission, 
                             sojourn_shape = sc_hsmm_fit$model$sojourn$shape,
                             sojourn_scale = sc_hsmm_fit$model$sojourn$scale)))
    
  }
  } else {
    exp_hmm_cell.function <- function(single_cell){
      
      sc <- single_cell[,.(chr, start = fifelse(strand == 1, bp, bp - readlength), 
                           end   = fifelse(strand == 1, bp + readlength, bp)) 
      ][order(chr, start * 0.5 + end * 0.5)
      ][, dist := c(0, diff(start * 0.5 + end * 0.5)), .(chr)]
      
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
      
      list(data = data.table:::merge.data.table(sc, by = c("chr", "start", "end"), all.x = T,
                                                cbind(sc[dist > 0][,.(chr, start, end)], 
                                                      data.table(state     = sc_hsmm_predict$s, 
                                                                 posterior = sc_hsmm_predict$p[,2]))),
           model = unlist(list(init = sc_hsmm_fit$model$init, 
                               emission = sc_hsmm_fit$model$parms.emission, 
                               sojourn_shape = sc_hsmm_fit$model$sojourn$shape,
                               sojourn_scale = sc_hsmm_fit$model$sojourn$scale)))
      
    }
  }
  
  
  {if(is_allelic){
    fread(list.files(pattern = "_filterer.tsv.gz"))[,.(cell, chr, bp, allele, readlength, strand)]
  } else {
    fread(list.files(pattern = "_filterer.tsv.gz"))[,.(cell, chr, bp, readlength, strand)]
  }} %>%
    data.table:::split.data.table(by = "cell", keep.by = F) %>%
    map(exp_hmm_cell.function) %>% #, mc.cores = cpus) %>% 
    {
      rbindlist(map(.[map_lgl(., is.list)], `[[`, "data"), idcol = "cell") %>%
        write_tsv(file = paste0(exp, "_hmmfitdata.tsv.gz"))
      rbindlist(map(.[map_lgl(., is.list)], function(x) data.table(parameter = names(x$model), value = x$model)), idcol = "cell") %>%
        write_tsv(file = paste0(exp, "_hmmfitmodel.tsv.gz"))
    }
  
  
} #number 2

super_singlefork   <- function(exp, is_allelic = TRUE){
  
  setDTthreads(threads = 0, restore_after_fork = T)
  
  try({for(i in 1:nrow(config)) assign(config$V1[i], config$V2[i])}, silent = T)
  
  if(is_allelic) {
  
    fread(list.files(pattern = "_hmmfitdata.tsv.gz")
    )[order(cell, chr, start * 0.5 + end * 0.5)
    ][, state_after := lead(state), .(cell, chr)
    ][is.na(state), state := state_after
    ][is.na(state_after), state_after := state  
    ][, group := cumsum((state == 1 & state_after == 2) | 
                          (state == 2 & state_after == 1)), .(cell, chr)
    ][state == 2 & state_after == 1, group := group - 1
    ][,.(
      state = max(state),
      allele_asigned = fcase(("A" %in% allele & "B" %in% allele) | ("both" %in% allele), "both",
                             "A" %in% allele, "A",
                             "B" %in% allele, "B", default = "none"),
      N = .N,
      center = min(start) * 0.5 + max(end) * 0.5,
      width = max(end) - min(start) + 1,
      coverage = sum(end - start)),
      .(cell, chr, group)] %>%
      write_tsv(file = paste0(exp, "_singlefork.tsv.gz"))
  } else {
    
    fread(list.files(pattern = "_hmmfitdata.tsv.gz")
    )[order(cell, chr, start * 0.5 + end * 0.5)
    ][, state_after := lead(state), .(cell, chr)
    ][is.na(state), state := state_after
    ][is.na(state_after), state_after := state  
    ][, group := cumsum((state == 1 & state_after == 2) | 
                          (state == 2 & state_after == 1)), .(cell, chr)
    ][state == 2 & state_after == 1, group := group - 1
    ][,.(
      state = max(state),
      allele_asigned = "none",
      N = .N,
      center = min(start) * 0.5 + max(end) * 0.5,
      width = max(end) - min(start) + 1,
      coverage = sum(end - start)),
      .(cell, chr, group)] %>%
      write_tsv(file = paste0(exp, "_singlefork.tsv.gz"))
    
  }
  
} #number 2

# actually done these jobs

all_jobs <- c()


##### job submitting part

# 1 poissontest
cat("\n--\n----- 1.poissontest -----\n--\n")
if(!file_check["poissontest"]){
  
  poissontest <- Slurm_EvalQ_custom(job_name = paste0(jobid, ".poissontest"), njobs = 1, plan = "wait", 
                                    expr = {super_poissontest(jobid)}, 
                                    sbatch_opt = list(mem = "10G", `cpus-per-task` = "4", time = "1:00:00"), 
                                    export = c("super_poissontest", "jobid", "config"))
  
  #print(readRDS(list.files(paste0(jobid, ".poissontest/"), pattern = "03", full.names = T)))
  
  Slurm_clean(poissontest)
  
  file_check <- file.checker()
  
  all_jobs["poissontest"] <- poissontest$jobid
  
}


# 2 filterer
cat("\n--\n----- 1.filterer -----\n--\n")
if(file_check["poissontest"] & # dependency
   !file_check["filterer"]){
  filterer <- Slurm_EvalQ_custom(job_name = paste0(jobid, ".filterer"), njobs = 1, plan = "wait", 
                                 expr = {super_filterer(jobid)}, 
                                 sbatch_opt = list(mem = "20G", `cpus-per-task` = "4", time = "2:00:00"), 
                                 export = c("super_filterer", "jobid", "config"))
  
  print(readRDS(list.files(paste0(jobid, ".filterer/"), pattern = "03", full.names = T)))
  
  Slurm_clean(filterer)
  
  file_check <- file.checker()
  
  all_jobs["filterer"] <-filterer$jobid
  
}


# 3 overlapscore
cat("\n--\n----- 3.overlapscore -----\n--\n")
if(file_check["filterer"] &  # dependency
   !file_check["overlapscore"]){
  
  input_size <- file.size(list.files(pattern = "_filterer")) / 1e9
  
  for(attempt in 1:3){
    
    est_mem <- paste0(min(max(c(round(input_size * 100 * attempt), 10)), 200), "G")
    
    print(paste(attempt, est_mem))
    
    overlapscore <- Slurm_EvalQ_custom(job_name = paste0(jobid, ".overlapscore"), njobs = 1, plan = "wait", 
                                       expr = {super_overlapscore(jobid)}, 
                                       sbatch_opt = list(mem = est_mem, `cpus-per-task` = "4", time = "10:00:00"), 
                                       export = c("super_overlapscore", "jobid", "config"))
    
    #print(readRDS(list.files(paste0(jobid, ".overlapscore/"), pattern = "03", full.names = T)))
    
    Slurm_clean(overlapscore)
    
    if(as.numeric(status(overlapscore)) == 0) break
    
  }
  
  file_check <- file.checker()
  
  all_jobs["overlapscore"] <- overlapscore$jobid
  
}


# 4 sphaseorder
cat("\n--\n----- 4.sphaseorder -----\n--\n")
if(file_check["overlapscore"] &  # dependency
   !file_check["sphaseorder"]){
  
  sphaseorder <- Slurm_EvalQ_custom(job_name = paste0(jobid, ".sphaseorder"), njobs = 1, plan = "wait", 
                                    expr = {super_sphaseorder(jobid)}, 
                                    sbatch_opt = list(mem = "10G", `cpus-per-task` = "2", time = "2:00:00"), 
                                    export = c("super_sphaseorder", "jobid", "config"))
  
  #print(readRDS(list.files(paste0(jobid, ".sphaseorder/"), pattern = "03", full.names = T)))
  
  Slurm_clean(sphaseorder)
  
  file_check <- file.checker()
  
  all_jobs["sphaseorder"] <- sphaseorder$jobid
  
}


# 5 pc
cat("\n--\n----- 5.pc -----\n--\n")
if(file_check["filterer"] &  # dependency
   !file_check["pc"]){
  
  input_size <- file.size(list.files(pattern = "_filterer")) / 1e9
  
  for(attempt in 1:3){
    
    est_mem <- paste0(min(max(c(round(input_size * 75 * attempt), 10)), 200), "G")
    
    print(paste(attempt, est_mem))
    
    pc <- Slurm_EvalQ_custom(job_name = paste0(jobid, ".pc"), njobs = 1, plan = "wait", 
                             expr = {super_pc(jobid)}, 
                             sbatch_opt = list(mem = est_mem, `cpus-per-task` = "4", time = "2:00:00"), 
                             export = c("super_pc", "jobid", "path2source", "config"))
    
    #print(readRDS(list.files(paste0(jobid, ".pc/"), pattern = "03", full.names = T)))
    
    Slurm_clean(pc)
    
    if(as.numeric(status(pc)) == 0) break
    
  }
  
  file_check <- file.checker()
  
  all_jobs["pc"] <- pc$jobid
  
}


# 6 psmashc
cat("\n--\n----- 6.psmashc -----\n--\n")
if(file_check["pc"] &  # dependency
   !file_check["psmashc"]){
  
  input_size <- file.size(list.files(pattern = "_pc")) / 1e9
  
  for(attempt in 1:3){
    
    est_mem <- paste0(min(max(c(round(input_size * 7.5 * attempt), 10)), 200), "G")
    
    print(paste(attempt, est_mem))
    
    psmashc <- Slurm_EvalQ_custom(job_name = paste0(jobid, ".psmashc"), njobs = 1, plan = "wait", 
                             expr = {super_psmashc(jobid)}, 
                             sbatch_opt = list(mem = est_mem, `cpus-per-task` = "4", time = "2:00:00"), 
                             export = c("super_psmashc", "jobid", "config"))
    
    #print(readRDS(list.files(paste0(jobid, ".pc/"), pattern = "03", full.names = T)))
    
    Slurm_clean(psmashc)
    
    if(as.numeric(status(psmashc)) == 0) break
    
  }
  
  file_check <- file.checker()
  
  all_jobs["psmashc"] <- psmashc$jobid
  
}


# 7 mmfit
cat("\n--\n----- 7.mmfit -----\n--\n")
if(file_check["pc"] &  # dependency
   !file_check["mmfit"]){
  
  input_size <- file.size(list.files(pattern = "_pc")) / 1e9
  
  for(attempt in 1:3){
    
    est_mem <- paste0(min(max(c(round(input_size * 22 * attempt), 10)), 200), "G")
    
    print(paste(attempt, est_mem))
    
    mmfit <- Slurm_EvalQ_custom(job_name = paste0(jobid, ".mmfit"), njobs = 1, plan = "wait", 
                                expr = {super_mmfit(jobid)}, 
                                sbatch_opt = list(mem = est_mem, `cpus-per-task` = "4", time = "15:00:00"), 
                                export = c("super_mmfit", "jobid", "path2source", "config"))
    
    #print(readRDS(list.files(paste0(jobid, ".mmfit/"), pattern = "03", full.names = T)))
    
    Slurm_clean(mmfit)
    
    if(as.numeric(status(mmfit)) == 0) break
    
  }
  
  file_check <- file.checker()
  
  all_jobs["mmfit"] <- mmfit$jobid
  
}


# 8 BICfit
cat("\n--\n----- 8.BICfit -----\n--\n")
if(file_check["pc"] & file_check["mmfit"] & # dependency
   !file_check["BICfit"]){
  
  input_size <- file.size(list.files(pattern = "_pc")) / 1e9
  
  for(attempt in 1:3){
    
    est_mem <- paste0(min(max(c(round(input_size * 10 * attempt), 10)), 200), "G")
    
    print(paste(attempt, est_mem))
  
  BICfit <- Slurm_EvalQ_custom(job_name = paste0(jobid, ".BICfit"), njobs = 1, plan = "wait", 
                               expr = {super_BICfit(jobid)}, 
                               sbatch_opt = list(mem = est_mem, `cpus-per-task` = "2", time = "1:00:00"), 
                               export = c("super_BICfit", "jobid", "path2source", "config"))
  
  #print(readRDS(list.files(paste0(jobid, ".BICfit/"), pattern = "03", full.names = T)))
  
  Slurm_clean(BICfit)
  
  if(as.numeric(status(BICfit)) == 0) break
  
  }
  
  file_check <- file.checker()
  
  all_jobs["BICfit"] <- BICfit$jobid
  
}


# 9 hmmfit
cat("\n--\n----- 9.hmmfit -----\n--\n")
if(file_check["filterer"] &  # dependency
   !file_check["hmmfit"]){
  
  input_size <- file.size(list.files(pattern = "_filterer")) / 1e9
  
  for(attempt in 1:3){
    
    est_mem <- paste0(min(max(c(round(input_size * 30 * attempt), 10)), 200), "G")
    
    print(paste(attempt, est_mem))
    
    hmmfit <- Slurm_EvalQ_custom(job_name = paste0(jobid, ".hmmfit"), njobs = 1, plan = "wait", 
                                 expr = {super_hmmfit(jobid)}, 
                                 sbatch_opt = list(mem = est_mem, `cpus-per-task` = "4", time = "10:00:00", gres = "tmpspace:100G"), 
                                 export = c("super_hmmfit", "jobid", "config"))
    
    print(readRDS(list.files(paste0(jobid, ".hmmfit/"), pattern = "03", full.names = T)))
    
    Slurm_clean(hmmfit)
    
    if(as.numeric(status(hmmfit)) == 0) break
    
  }
  
  file_check <- file.checker()
  
  all_jobs["hmmfit"] <- hmmfit$jobid
  
}


# 10 singlefork
cat("\n--\n----- 10.singlefork -----\n--\n")
if(file_check["hmmfit"] &  # dependency
   !file_check["singlefork"]){
  
  singlefork <- Slurm_EvalQ_custom(job_name = paste0(jobid, ".singlefork"), njobs = 1, plan = "wait", 
                                   expr = {super_singlefork(jobid)}, 
                                   sbatch_opt = list(mem = "10G", `cpus-per-task` = "4", time = "2:00:00"), 
                                   export = c("super_singlefork", "jobid", "config"))
  
  print(readRDS(list.files(paste0(jobid, ".singlefork/"), pattern = "03", full.names = T)))
  
  Slurm_clean(singlefork)
  
  file_check <- file.checker()
  
  all_jobs["singlefork"] <- singlefork$jobid
  
}



system(paste0("sacct -o JobID,JobIDRaw,JobName,AllocCPUS,ReqMem,MaxRSS,timelimit,Elapsed,State -j ",
              paste0(all_jobs, collapse = ","),
              " -P > ", jobid, "_", gsub("-| |:", "", Sys.time()), "_memstats.tsv"))

warnings()


cat("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     \n END END END END END END END END END END 
     \n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")





