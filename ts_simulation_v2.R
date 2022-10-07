library(slendr)
library(tidyverse)
library(sf)
library(cluster)
library(factoextra)
library(lattice)

model_dir <- tempdir()
source("~/project/scripts/parameters_grid_function.R")

grid <- parameters_grid(5, c(5, 7.5, 10), 10, 10)
grid <- cbind(grid, percent_not_coal = 0.0, recap_diversity = NA)
grid

ts_dir <- "ts_out5"
# in put1 ci sono i recapitated ts, 400 den 5, r 5
#same in output2 ma con le diversity giuste

dir.create(ts_dir)


for (i in seq_len(nrow(grid))) {
  radius <-   grid[i, ]$radius
  N <-        grid[i, ]$N
  density <-  grid[i, ]$density
  competition <-  grid[i, ]$competition
  dispersal <- grid[i,  ]$dispersal
  mating <-   grid[i, ] $mating
  burnin <-   grid[i, ]$burnin
  print(paste0("simulating:   ", i))
  
  map <- world(xrange = c(0, 2 * radius), 
          yrange = c(0, 2 * radius), landscape = "blank")

  pop <- population("pop", time = 1, N = N, map = map,
          center = c(radius, radius), radius = radius)

  ts_file <- file.path(ts_dir, sprintf("density%s_comp%s_disp%s.trees",
       density, competition, dispersal))
  
  model <- compile_model(populations = pop, resolution = 1, generation_time = 1,
          simulation_length = 5000, competition = competition,
          dispersal = dispersal, mating = mating,
          path = model_dir, overwrite = TRUE, force = TRUE)
  
  possible_error <- tryCatch(
          ts <- slim(model, sequence_length = 50e6,
          recombination_rate = 1e-8, random_seed = 1810,
          #samples = t,
          burnin = burnin * N, output = ts_file, load = TRUE),
  error = function(e) e
  )
  if (inherits(possible_error, "error")) next
  
  ts <- ts_mutate(ts, mutation_rate = 1e-6)
  sample <- ts_samples(ts)
  sets <- split(sample, sample$time) %>% lapply(function(pop) pop$name)
  overall_diversity_over_time <- ts_diversity(ts, sets)
  grid[i, ]$pop_diversity <- mean(overall_diversity_over_time$diversity)
  
  if (ts_coalesced(ts) == "FALSE") {
        grid[i, ]$coalesced <- FALSE
        perc <- length(ts_coalesced(ts, TRUE)) / ts$num_trees
        grid[i, ]$percent_not_coal <- perc * 100 }
  
  
  ts <- ts_recapitate(ts, recombination_rate = 1e-8, Ne = N)
  grid[i, ]$recap_diversity <- mean(ts_diversity(ts, sets)$diversity)


  }

grid

write_tsv(grid, "~/project/Thesis/den7_5.tsv")




