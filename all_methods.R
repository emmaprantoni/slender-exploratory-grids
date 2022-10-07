library(slendr)
library(tidyverse)
library(sf)
library(cluster)
library(factoextra)
library(lattice)
library(viridis)


grid <- read_tsv("~/project/Thesis/analysed_recp_400.tsv")
ts_dir <- "ts_out_den5_b10x_400"

grid


grid <- cbind(grid,
        kmeans = NA,
        hclust = NA,
        agnes = NA,
        diana = NA
        )

for (i in seq_len(nrow(grid))){
    cat(sprintf("analyzing simulation from grid %d/%d\n", i, nrow(grid)))
    den <- grid[i, ]$density
    competition <-  grid[i, ]$competition
    dispersal <-  grid[i, ]$dispersal
    mating <-  grid[i, ]$mating

    ts_file <- file.path(ts_dir, sprintf("density%s_comp%s_disp%s.trees",
        density, competition, dispersal))

    possible_error <- tryCatch(
        ts <- ts_load(ts_file, simplify = TRUE,
        mutate = TRUE, mutation_rate = 1e-6),
        error = function(e) e)
    if (inherits(possible_error, "error")) next
 
    nodes <- ts_nodes(ts) %>% 
    filter(time == min(time)) %>%
    distinct(ind_id, .keep_all = TRUE)

    points <- st_coordinates(nodes, scale = TRUE)
    gs <- clusGap(points, K.max = 15, FUN = kmeans, nstart = 50, B = 50)
    k <-  maxSE(f         = gs$Tab[, "gap"],
            SE.f      = gs$Tab[, "SE.sim"],
            SE.factor = 1)
    grid[i, ]$kmeans <- k
    
    for (method in c('hclust', 'agnes', 'diana')){
    gs <- clusGap(points, K.max = 10,FUN = hcut, hc_func = method, hc_method = 'average', B = 20)
    k <-  maxSE(f         = gs$Tab[,"gap"],
                SE.f      = gs$Tab[,"SE.sim"],
                SE.factor = 1)
    grid[i, method] <- k
  }
}

grid

pdf('~/project/methods_plots.pdf')
for ( cl in list (kmeans, hclust, agnes, diana)) {
p1 <- levelplot(kmeans ~ competition * dispersal, grid, col.regions = rainbow(16),
        at = 0:16, main = str(cl))
print(p1) }
dev.off()


