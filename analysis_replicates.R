grid <- read_tsv("~/project/Thesis/den_5x5_b10x_rep5_full.tsv")
ts_dir <- "ts_out_den5_5x5_b20x"

grid


grid <- cbind(grid,
        km_clusters = NA,
        hc_clusters = NA,
        km_within_div = NA,
        hc_within_div = NA,
        km_median_size = NA,
        hc_median_size = NA,
        km_Ne = NA,
        hc_Ne = NA,
        pop_Ne = NA,
        km_mode = NA,
        hc_mode = NA)



for (i in seq_len(nrow(grid))){
    km_c <- rep(NA,5)
    hc_c <- rep(NA,5)
    km_w <- rep(NA,5)
    hc_w <- rep(NA,5)
    km_s <- rep(NA,5)
    hc_s <- rep(NA,5)
    km_n <- rep(NA,5)
    hc_n <- rep(NA,5)
    pop_n <- rep(NA,5)

    cat(sprintf("analyzing simulation from grid %d/%d\n", i, nrow(grid)))
    den <- grid[i, ]$density
    competition <-  grid[i, ]$competition
    dispersal <-  grid[i, ]$dispersal
    mating <-  grid[i, ]$mating

    for(j in 1:5) {
        ts_file <- file.path(ts_dir, sprintf("density%s_comp%s_disp%s_rep%s.trees",
        density, competition, dispersal,j))


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
    km_c[j] <- k

    km <- hcut(points, k, hc_func = "hclust", hc_method = "average")
    #perchè il clustering funziona meglio
    km_s[j] <- median(table(km$cl))

    km_tib <- tibble(cbind(nodes, km_assigment = km$cluster))
    sample_set <- split(km_tib, km_tib$km_assigment) %>%
                    lapply(function(pop) pop$node_id)
    result <- ts_diversity(ts, sample_sets = sample_set)
    km_w[j] <- mean(result$diversity)



    gs <- clusGap(points, K.max = 15, FUN = hcut, hc_func = "hclust",
        hc_method = "average",  B = 50)

    k <-  maxSE(f         = gs$Tab[, "gap"],
            SE.f      = gs$Tab[, "SE.sim"],
            SE.factor = 1)

    hc_c[j] <- k
    hc <- hcut(points, k, hc_func = "hclust", hc_method = "average")
    hc_s[j]<- median(table(hc$cl))

    hc_tib <- tibble(cbind(nodes, hc_assigment = hc$cluster))
    sample_set <- split(hc_tib, hc_tib$hc_assigment) %>%
                    lapply(function(pop) pop$node_id)
    result <- ts_diversity(ts, sample_sets = sample_set)

    hc_w[j] <- mean(result$diversity)
    
    }

    grid[i,]$km_clusters <- mean(km_c)
    grid[i,]$hc_clusters <- mean(hc_c)
    grid[i,]$km_within_div <- mean(km_w)
    grid[i,]$hc_within_div <- mean(hc_w)
    grid[i,]$km_median_size <- mean(km_s)
    grid[i,]$hc_median_size <- mean(hc_s)
    grid[i,]$km_Ne <- round(grid[i, ]$km_within_div / (4 * 1e-6))
    grid[i,]$hc_Ne <- round(grid[i, ]$hc_within_div / (4 * 1e-6))
    grid[i,]$pop_Ne <- round(grid[i, ]$pop_diversity / (4 * 1e-6))





}



grid


write_tsv(grid, "~/project/Thesis/analysed_replicates_recapitation_400.tsv")
