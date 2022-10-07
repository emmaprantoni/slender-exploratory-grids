library(slendr)
library(tidyverse)
library(sf)
library(cluster)
library(factoextra)
library(lattice)

model_dir <- tempdir()
source("~/project/scripts/parameters_grid_function.R")

grid <- parameters_grid(5, c(5), 10, 20)
grid <- cbind(grid, percent_not_coal = 0.0, recap_diversity = NA)
grid

ts_dir <- "replicates_recapitation_400"

dir.create(ts_dir)


for (i in seq_len(nrow(grid))) {
  radius <- grid[i,]$radius
  N <- grid[i,]$N
  density <- grid[i,]$density
  competition <- grid[i,]$competition
  dispersal <- grid[i,]$dispersal
  mating <- grid[i,]$mating
  burnin <- grid[i, ]$burnin
  print(paste0('simulating:  ', i))
  
  map <- world(xrange = c(0,2*radius), yrange = c(0,2*radius), landscape = 'blank')
  pop <- population('pop', time = 1, N = N, map = map, center = c(radius,radius), 
                    radius = radius)


  
  
  model <- compile_model(populations = pop, resolution = 1, generation_time = 1,
                         simulation_length = 5000, competition = competition, 
                         dispersal = dispersal, mating = mating,
                         path = model_dir, overwrite = TRUE, force = TRUE)
  
  #t <- schedule_sampling(model, times = c(), list(pop, N))
  
   rep_percent <- rep(NA,5)
   rep_diversity <- rep(NA,5)
   recap_div <- rep(NA,5)

  for(j in 1:5) {
        ts_file <- file.path(ts_dir, sprintf("density%s_comp%s_disp%s_rep%s.trees",
        density, competition, dispersal,j))
        possible_error <- tryCatch(
    ts <- slim(model, sequence_length = 50e6, recombination_rate = 1e-8, random_seed = 1810,
               #samples = t, 
               burnin = burnin*N, output = ts_file, load = TRUE),
    error = function(e) e
  )
  if (inherits(possible_error, "error")) next
  
  ts <- ts_mutate(ts, mutation_rate = 1e-6)
  sample <- ts_samples(ts)
  sets <- split(sample, sample$time)  %>% lapply(function(pop) pop$name)
  overall_diversity_over_time <- ts_diversity(ts, sets)
  rep_diversity[j] <- overall_diversity_over_time$diversity
  if(ts_coalesced(ts) == 'FALSE') grid[i,]$coalesced <- FALSE
  perc <- (length(ts_coalesced(ts, T))/ts$num_trees)*100
  rep_percent[j] <- perc

  if (ts_coalesced(ts) == "FALSE") {
        grid[i, ]$coalesced <- FALSE
        perc <- (length(ts_coalesced(ts, TRUE)) / ts$num_trees)*100
        rep_percent[j] <- perc }
  
  
  ts <- ts_recapitate(ts, recombination_rate = 1e-8, Ne = N) 
  recap_div[j] <- mean(ts_diversity(ts, sets)$diversity)




  }
  
  grid[i, ]$pop_diversity <- mean(rep_diversity)
  grid[i, ]$percent_not_coal <- mean(rep_percent)
  grid[i, ]$recap_diversity <- mean(recap_div)
  }


grid