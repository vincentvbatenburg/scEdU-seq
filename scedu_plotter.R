#!/usr/bin/env Rscript

# dependency

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

mixdist   <- function(Estimate, diffusion){
  
  if(!diffusion){
    if(length(Estimate) == 4){
      
      X <- data.table(dist = seq(0, 4e5, 5000)
      )[,.(dist = dist,
           E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[3]},
           U = dunif(dist, 0, 4e5) %>% {. / sum(.) * Estimate[4]} )] %>%
        melt(id.vars = "dist")
    } else
      if(length(Estimate) == 7){
        
        X <-  data.table(dist = seq(0, 4e5, 5000)
        )[,.(dist = dist,
             E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[5]} ,
             M = dnorm(dist, Estimate[2], sqrt(Estimate[3])) %>% {. / sum(.) * Estimate[6] },
             U = dunif(dist, 0, 4e5) %>% {. / sum(.) * Estimate[7]})] %>%
          melt(id.vars = "dist")
      } else 
        if(length(Estimate) == 10){
          X <-   data.table(dist = seq(0, 4e5, 5000)
          )[,.(dist = dist,
               E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[7]} ,
               M0 = dnorm(dist, Estimate[2], sqrt(Estimate[3])) %>% {. / sum(.) * Estimate[8] },
               M = dnorm(dist, Estimate[4], sqrt(Estimate[5])) %>% {. / sum(.) * Estimate[9] },
               U = dunif(dist, 0, 4e5) %>% {. / sum(.) * Estimate[10]})] %>%
            melt(id.vars = "dist")}
  } 
  
  if(diffusion){
    if(length(Estimate) == 4){
      
       X <- data.table(dist = seq(0, 4e5, 5000)
      )[,.(dist = dist,
           E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[3]},
           U = dunif(dist, 0, 4e5) %>% {. / sum(.) * Estimate[4]} )] %>%
        melt(id.vars = "dist")
      
      
    } else 
      if(length(Estimate) == 7){
        
        X <- data.table(dist = seq(0, 10e5, 5000)
        )[,.(dist = dist,
             E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[5]} ,
             M0 = dnorm(dist, Estimate[2], sqrt(Estimate[3])) %>% {. / sum(.) * Estimate[6] },
             U = dunif(dist, 0, 10e5) %>% {. / sum(.) * Estimate[7]})] %>%
          melt(id.vars = "dist")
        
      } else 
        if(length(Estimate) == 10){
          X <-   data.table(dist = seq(0, 10e5, 5000)
          )[,.(dist = dist,
               E = dexp(dist, 1 / Estimate[1]) %>% {. / sum(.) * Estimate[7]} ,
               M0 = dnorm(dist, Estimate[2], sqrt(Estimate[3])) %>% {. / sum(.) * Estimate[8] },
               M = dnorm(dist, Estimate[4], sqrt(Estimate[5])) %>% {. / sum(.) * Estimate[9] },
               U = dunif(dist, 0, 10e5) %>% {. / m(.) * Estimate[10]})] %>%
            melt(id.vars = "dist")}
  } 
  
  
  return(X[,variable := factor(variable, c("M", "E","M0", "U"))][])
  
  
  
}
frollconf <- function(X, Y, N){
  Ypad <- c(rep(NA, round(1.5 * N)), Y, rep(NA, round(1.5 * N)))
  
  data.table(ymeanpad = frollmean(Ypad, n = 2 * N + 1, align = "center", na.rm = T, hasNA = T)[],
             y95pad = frollapply(Ypad, n = 2 * N + 1, FUN = sd, align = "center", na.rm = T) * 2
  )[,.(x = X,
       ymean = ymeanpad[!is.na(Ypad)],
       y95 = y95pad[!is.na(Ypad)],
       yse = c(y95pad / sqrt(frollsum(!is.na(Ypad) * 1, n = 2 * N + 1, align = "center")))[!is.na(Ypad)])
  ]
  
}

dependencies <- setNames(c("_filterer",    "_overlapscore", "_filterer", "_pc",      "_pc",    "_pc", "_filterer", "_hmmfit"),
                         c("overlapscore", "sphaseorder", "pc",       "psmashc", "mmfit", "BICfit", "hmmfit", "singlefork"))

# arguments

arguments <- commandArgs(trailingOnly=TRUE)

print(arguments)

exp_pattern <- arguments[1] #"CD8"

base_folder <- paste0(getwd(), "/")

#if(length(arguments) > 1) base_folder <- arguments[2] #"hpc2/vincentvb/tcell_parp/EdU/"

sphase_flip <- "emptyplaceholder" #"DMSO-2hr|V1-2hr"

if(length(arguments) > 1) sphase_flip <- arguments[2] #"hpc2/vincentvb/tcell_parp/EdU/"

#
all_plots <- list()

for(i in 1:26) all_plots[[paste0(LETTERS, LETTERS)[i]]] <- "empty"

##### PLOT files
cat("\n1. PLOT files \n")
try({all_plots$AA <- data.table(file = list.files(base_folder, pattern = exp_pattern, recursive = T, include.dirs = F)
)[, exp := gsub(".*/|_.*|-pl.*","", file)
][, type := sub("_|-", "", gsub(paste0(exp, "|.tsv.gz|.RDS"), "", sub(".*/", "", file))), .(file)
][, size := file.size(paste0(base_folder, file)) / 1e6
][str_detect(type, "memstats|plotted", T) & str_detect(exp, "plot", T)] %>%
  ggplot() + geom_raster(aes(x = type, y = exp, fill = (size))) +
  facet_grid(cols = vars(str_detect(type, "pl")), scales = "free")
})


##### PLOT memory usage
cat("\n2. PLOT memory usage \n")
try({all_plots$BB <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*memstats"), recursive = T), 
                         function(x){fread(paste0(base_folder, x))[, exp := gsub(".*/|_.*", "", x)]})
)[, job := sub("_.*", "", JobID)
][, jobtype := sub(".*\\.", "", JobName)
][,.(mem_used = sum(as.numeric(sub("K", "", MaxRSS)), na.rm = T) / 1e6,
     time_req = Timelimit[1], 
     jobtype = jobtype[1],
     time_elapsed = Elapsed[1]), .(exp, job, AllocCPUS, mem_req = as.numeric(sub("Gn", "", ReqMem)))
][, dependency := grep(dependencies[jobtype], 
                       list.files(base_folder, pattern = exp, full.names = T, recursive = T), value = T)[1], .(exp, job)
][, mem_depend := file.size(dependency) / 1e9
][!is.na(mem_depend), depend_ratio := coef(lm(mem_used ~ mem_depend, data = .SD))[2], .(jobtype)
][] %>%
  ggplot() + 
  #geom_point(aes(x = as.POSIXct(time_req), y = as.POSIXct(time_elapsed), col = jobtype)) + 
  geom_point(aes(x = mem_depend, y = mem_used, col = jobtype)) + 
  geom_smooth(aes(x = mem_depend, y = mem_used, col = jobtype), method = "lm", se = F) + 
  geom_label(aes(x = 0.6, y = as.numeric(as.factor(jobtype)) * 2, label = paste(jobtype, round(depend_ratio, 2)), col = jobtype)) +
  #geom_point(aes(x = mem_req, y = mem_used, col = jobtype)) + 
  geom_abline(slope = 1, intercept = 0)})


##### PLOT poisson test
cat("\n3. PLOT poisson test \n")
try({all_plots$CC <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*poissontest"), recursive = T), 
                         function(x){fread(paste0(base_folder, x))[, exp := gsub(".*/|_.*", "", x)]})
) %>%
  ggplot() +
  geom_point(aes(x = log(mean), y = poisson_score, col = exp),
             size = 0.1) +
  facet_wrap(vars(exp), scales = "free") +
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0.2) +
  geom_vline(xintercept = c(-2.5, 2.5))})


##### PLOT posterior per experiment
cat("\n4. PLOT posterior per experiment \n")
try({all_plots$DD <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*filterer"), recursive = T), 
                         function(x){fread(paste0(base_folder, x))[, exp := gsub(".*/|_.*", "", x)]})
)[,.(N = as.double(.N)), .(exp, posterior)
][,N := N / max(N), .(exp)] %>%
  ggplot() + geom_tile(aes(x = 1, y = posterior, height = 1, width = N)) +
  facet_wrap(vars(exp))})


##### PLOT shpase order QCs
cat("\n5. PLOT shpase order QCs \n")
try({all_plots$EE <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*sphaseorder"), recursive = T), 
                         function(x){readRDS(paste0(base_folder, x))$bootstrap[, exp := gsub(".*/|_.*", "", x)]})
) %>% 
  ggplot() + geom_point(aes(x = rank, y = value, col = N), size = 0.2) + 
  facet_wrap(vars(exp)) +
  scale_color_viridis_c()
})

try({all_plots$FF <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*sphaseorder"), recursive = T), 
                         function(x){as.data.table(readRDS(paste0(base_folder, x))$final_order)[, exp := gsub(".*/|_.*", "", x)]})
) %>% ggplot(aes(y = X1.1, x = X2.1)) +
  geom_segment(aes(yend = X2, xend = Y2), size = 0.05) +
  geom_point(aes(col = rank), size = 1) +
  coord_fixed() +
  facet_wrap(vars(exp)) +
  scale_color_viridis_c(name = "Tour index")
})

####### plotting cut tracks
cat("\n6. PLOT cut tracks \n")
try({all_plots$GG <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*filterer"), recursive = T), 
                         function(x){
                           fread(paste0(base_folder, x)
                           )[chr == 2 & posterior > 900
                           ][, exp := gsub(".*/|_.*", "", x)
                           ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                           )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                           ][, .(N = as.double(.N)), .(exp, rank, bp = round(bp, -5))
                           ][, N := (N) / quantile(N, 0.95), .(rank)
                           ][bp > 3e7 & bp < 9e7
                           ][N > 1, N := 1]
                         })
)[str_detect(exp, sphase_flip), rank := abs(rank - 1)  #"DMSO-2hr|V1-2hr"
]%>% 
  ggplot() + 
  geom_raster(aes(y = rank, x = bp / 1e6, fill = (N)) ) +
  scale_y_continuous(trans = "reverse") + ylab("S-phase Progression") + xlab("Chromosome 2 (Mbp)") +
  facet_wrap(vars(exp), scales = "free", ncol = 2) +
  scale_fill_gradientn(colours = c("white", "blue", "blue")) +
  theme_bw() +
  coord_cartesian(expand = F) })


####### plotting genome pile up
cat("\n7. PLOT genome pile up \n")
try({all_plots$HH <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*filterer"), recursive = T), 
                         function(x){
                           fread(paste0(base_folder, x)
                           )[chr %in% 1:22, .(N = as.double(.N), exp = gsub(".*/|_.*", "", x)), .(cell, chr, bp = round(bp / 1e6))
                           ][, N := N / sum(N), .(cell)
                           ][, .(N = sum(N)), .(exp, chr, bp)
                           ][, N := N / quantile(N, 0.97), .(chr)
                           ][N > 1, N := 1][]
                         })) %>% 
  ggplot() + 
  geom_col(aes(y = N, x = bp, fill = log10(N)) ) +
  facet_grid(rows = vars(exp), cols = vars(chr), scales = "free") +
  scale_fill_viridis_c() +
  coord_cartesian(expand = F)
})


####### plotting model fits
cat("\n8. PLOT model fits \n")
try({all_plots$II <- merge.data.table(
  rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                function(x){fread(paste0(base_folder, x)
                )[, exp := gsub(".*/|_.*", "", x)
                ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                ]})
  )[((!diffusion) & model != "quadruple") | (diffusion)
  ][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
  ][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
  ][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
  ][parameter != 11, mixdist(Estimate, diffusion = diffusion[1]), .(exp, cell, rank, diffusion)],
  rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*psmashc"), recursive = T), 
                function(x){fread(paste0(base_folder, x))[, exp := gsub(".*/|_.*", "", x)]})),
  by = c("exp", "cell", "dist")
)[dist > 0
][,  N := as.double(N)
][, N := N / sum(N), .(exp, cell, diffusion, variable)
][str_detect(exp, sphase_flip), rank := abs(rank - 1)   
  #][rank > 0.4 & rank < 0.6
][,.SD[cell %in% sample(unique(cell), 3)], .(exp)
][dist > 0 & dist < 4e5] %>%
  ggplot() + 
  geom_col(aes(x = dist, y = value, fill = variable)) + 
  geom_line(aes(x = dist , y = N)) +
  facet_wrap(vars(exp, rank), scales = "free", nrow = 4) 
})


##### PLOT speeds
cat("\n9. PLOT speeds \n")
try({all_plots$JJ <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                         function(x){fread(paste0(base_folder, x)
                         )[, exp := gsub(".*/|_.*", "", x)
                         ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                         )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                         ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp, sphase_flip), rank := abs(rank - 1)  
][rank < 0.8
  #][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% 
  {
    .[!is.na(rank)][order(rank), frollconf(rank, Estimate, 15), .(exp)]  %>% 
      ggplot() + 
      geom_point(data = ., aes(x = rank, y = Estimate, col = cell > 384)) + 
      geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95), alpha = 0.3) +
      geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse), alpha = 0.3) +
      geom_line(aes(x = x, y = ymean)) +
      facet_grid(cols = vars(exp)) + theme_bw() + ylab("Speed (kb/min)") + xlab("S-phase progression")
    
  }
})

try({all_plots$KK <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                         function(x){fread(paste0(base_folder, x)
                         )[, exp := gsub(".*/|_.*", "", x)
                         ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                         )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                         ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp, sphase_flip), rank := abs(rank - 1)  
][rank < 0.8
  
][!is.na(rank)][order(rank), frollconf(rank, Estimate, 15), .(exp)]  %>% 
  ggplot() + 
  #geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95, col = exp), alpha = 0.3) +
  geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse, fill = exp), alpha = 0.3) + 
  geom_line(aes(x = x, y = ymean, col = exp)) + 
  theme_bw() + ylab("Speed (kb/min)") + xlab("S-phase progression")
})

pc_sim <- map(1:210, function(x){ 
  if(length(list.files("/hpc/hub_oudenaarden/vincentvb/forkdot/pc_sim2/", paste0("_", x, ".tsv"))) > 0){
    fread(paste0("/hpc/hub_oudenaarden/vincentvb/forkdot/pc_sim2/pc_sim_super_", x, ".tsv"))
  }else{
    
    fread(paste0("/hpc/hub_oudenaarden/vincentvb/forkdot/pc_sim1/pc_sim_super_", x, ".tsv"))
  }
  
}) %>%
  rbindlist()

pc_sim[, c("u_R", "s_R")       := .(abs(u - mean(u)) / sd(u), abs(s - mean(s)) / sd(s)), .(U,S)
][, c("u_rank", "s_rank") := .(rank(u_R), rank(s_R)), .(U,S)]

U_model <- loess(U ~ u + s, pc_sim[u_rank < 10 | s_rank < 10])
S_model <- loess(S ~ u + s, pc_sim[u_rank < 10 | s_rank < 10])


##### PLOT corrected speeds
cat("\n10. PLOT speeds \n")
try({all_plots$LL <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                         function(x){fread(paste0(base_folder, x)
                         )[, exp := gsub(".*/|_.*", "", x)
                         ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                         )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                         ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp, sphase_flip), rank := abs(rank - 1)  
][rank < 0.8
][, u := predict(U_model, data.table(u = Estimate / 75000, s = est_sd / 75000))
] %>% 
  {
    .[!is.na(rank)][order(rank), frollconf(rank, u, 15), .(exp)]  %>% 
      ggplot() + 
      geom_point(data = ., aes(x = rank, y = u, col = cell > 384)) + 
      geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95), alpha = 0.3) +
      geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse), alpha = 0.3) +
      geom_line(aes(x = x, y = ymean)) +
      facet_grid(cols = vars(exp)) + theme_bw() + ylab("Speed (kb/min)") + xlab("S-phase progression")
    
  }
})

try({all_plots$MM <- rbindlist(map(list.files(base_folder, pattern = paste0(exp_pattern, ".*BICfit"), recursive = T), 
                         function(x){fread(paste0(base_folder, x)
                         )[, exp := gsub(".*/|_.*", "", x)
                         ][, rank := as.data.table(readRDS(list.files(base_folder, pattern = paste0(exp, "_sphaseorder"), recursive = T, full.names = T)
                         )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                         ]})
)[((!diffusion) & model != "quadruple") | (diffusion)
][log2(est_sd / stde_mean) > 2 | is.na(stde_mean)
][, failed := Estimate[parameter == 4] < 4e5 & est_sd < 6e4, .(exp, cell, diffusion, model, peak_init)
][(failed), .SD[BIC == min(BIC)], .(exp, diffusion, cell)
][parameter == 4 
][str_detect(exp, sphase_flip), rank := abs(rank - 1)  
][rank < 0.8
][, u := predict(U_model, data.table(u = Estimate / 90000, s = est_sd / 90000))
][!is.na(rank)][order(rank), frollconf(rank, u, 15), .(exp)]  %>% 
  ggplot() + 
  #geom_ribbon(aes(x = x, ymin = ymean - y95, ymax = ymean + y95, col = exp), alpha = 0.3) +
  geom_ribbon(aes(x = x, ymin = ymean - yse, ymax = ymean + yse, fill = exp), alpha = 0.3) + 
  geom_line(aes(x = x, y = ymean, col = exp)) + 
  theme_bw() + ylab("Speed (kb/min)") + xlab("S-phase progression")
})


saveRDS(all_plots, file = paste0("plotted_", exp_pattern, ".RDS"))

warnings()

