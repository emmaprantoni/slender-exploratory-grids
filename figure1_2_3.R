library(slendr)
library(sf)
library(tidyverse)
library(cluster)
library(ggplot2)
library(factoextra) # clustering visualization
library(dendextend)

set.seed(1810)
density  = 10
r = 5
N <- round(r^2*pi*density)
grid <- expand_grid(
  radius = r, 
  density = density,
  N = N, 
  dispersal  =c(0.4,0.2,0.1),
  competition = c(3,5,7,9),
  mating = 0.1,#c(0.5,2,5,10),
  k = NA
)
island <- region( map = map ,center = c(r,r), radius = r)
map <- world(xrange = c(0,2*r), yrange = c(0,2*r), landscape = island)

#pop1
pop <- population("pop1", time = 1, N = N, map = map, competition = 3)

model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 0.1 , dispersal = 0.1)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])
gs1 <- fviz_gap_stat(gap_stat, linecolor = rainbow(16)[6])
p1 <- ggplot() + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p1)
plot(gs1)
#pop2
pop <- population("pop1", time = 1, N = N, map = map, competition = 5)
model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 0.1 , dispersal = 0.1)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])

p2 <- ggplot() + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p2)


#pop3
pop <- population("pop1", time = 1, N = N, map = map, competition = 7)
model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 0.1 , dispersal = 0.1)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])

p3 <- ggplot() + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p3)

#pop4
pop <- population("pop1", time = 1, N = N, map = map, competition = 9)
model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 0.1 , dispersal = 0.1)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])

p4 <- ggplot() + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p4)




#pop1
pop <- population("pop1", time = 1, N = N, map = map, competition = 3)
model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 0.1 , dispersal = 0.2)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])
gs5 <- fviz_gap_stat(gap_stat, linecolor = rainbow(16)[4])
p5 <- ggplot() + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p5)

#pop2
pop <- population("pop1", time = 1, N = N, map = map, competition = 5)
model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 0.1 , dispersal = 0.2)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])

p6 <- ggplot() + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p6)


#pop3
pop <- population("pop1", time = 1, N = N, map = map, competition = 7)
model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 0.1 , dispersal = 0.2)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])
gs7 <- fviz_gap_stat(gap_stat, linecolor = 'darkorange1')
p7 <- ggplot() + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p7)

#pop4
pop <- population("pop1", time = 1, N = N, map = map, competition = 9)

model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 0.1 , dispersal = 0.2)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])
gs8 <- fviz_gap_stat(gap_stat, linecolor = 'brown1')
p8 <- ggplot() + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p8)





#pop1
pop <- population("pop1", time = 1, N = N, map = map, competition = 3)
model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 0.1 , dispersal = 0.4)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])

p9 <- ggplot() + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p9)

#pop2
pop <- population("pop1", time = 1, N = N, map = map, competition = 5)
model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 0.1 , dispersal = 0.4)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])

p10 <- ggplot() + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p10)


#pop3
pop <- population("pop1", time = 1, N = N, map = map, competition = 7)
model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 0.1 , dispersal = 0.4)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])

p11 <- ggplot() + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p11)

#pop4
pop <- population("pop1", time = 1, N = N, map = map, competition = 9)
model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 0.1 , dispersal = 0.4)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])

p12 <- ggplot(res = 1500) + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p12)



ggpubr::ggarrange(p9,p10,p11,p12,p5,p6,p7,p8,p1,p2,p3,p4,ncol= 4, nrow = 3)
ggpubr::ggarrange(p5,p6,p7,p8,p1,p2,p3,p4, ncol= 4, nrow = 2)
gs1 <- fviz_gap_stat(gap_stat, linecolor ='black')
plot(gs11)

ggpubr::ggarrange(p1,p5,p7,p8, ncol= 4)
ggpubr::ggarrange(gs1,gs5,gs7,gs8, ncol= 2, nrow=2)

pop <- population("pop1", time = 1, N = N, map = map, competition = 3)
model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 1 , dispersal = 0.1)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])

p13 <- ggplot() + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p13)

#pop2
pop <- population("pop1", time = 1, N = N, map = map, competition = 3)
model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 2 , dispersal = 0.1)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])

p14 <- ggplot() + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p14)


#pop3
pop <- population("pop1", time = 1, N = N, map = map, competition = 3)
model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 5 , dispersal = 0.1)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])

p15 <- ggplot() + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p15)

#pop4
pop <- population("pop1", time = 1, N = N, map = map, competition = 3)
model <- compile_model(pop, generation_time = 1, simulation_length = 5000,
                       resolution = 0.1, mating = 10 , dispersal = 0.1)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])

p16 <- ggplot(res = 1500) + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p16)


ggpubr::ggarrange(p13, p14,p15, p16, ncol= 4)





## Figure 3 ####
#uniformily distributed
pop <- population("pop1", time = 1, N = N, map = map, competition = 7)
model <- compile_model(pop, generation_time = 1, simulation_length = 1000,
                       resolution = 0.1, mating = 5 , dispersal = 5)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])
km <- hcut(points, k, hc_func = 'hclust', hc_method = 'average')
p1 <- ggplot(res = 1500) + geom_sf(data = map) + geom_sf(data = nodes, col = km$cluster+2)
plot(p1)
gs1 <- fviz_gap_stat(gap_stat, linecolor = 'cornflowerblue')



pop <- population("pop1", time = 1, N = N, map = map, competition = 4.5)
model <- compile_model(pop, generation_time = 1, simulation_length = 1000,
                       resolution = 0.1, mating = 0.5 , dispersal = 0.5)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])
k
p2 <- ggplot(res = 1500) + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p2)
gs2 <- fviz_gap_stat(gap_stat, linecolor = 'cornflowerblue')


pop <- population("pop1", time = 1, N = N, map = map, competition = 10)
model <- compile_model(pop, generation_time = 1, simulation_length = 1000,
                       resolution = 0.1, mating = 0.1 , dispersal = 0.1) #0.2

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])
hc <- hcut(points, 2, hc_func = 'hclust', hc_method = 'average')
p3 <- ggplot(res = 1500) + geom_sf(data = map) + geom_sf(data = nodes, col = ifelse(hc$cl == 1, 'palevioletred', 'gold2'))
plot(p3)
gs3 <- fviz_gap_stat(gap_stat, linecolor = 'palevioletred')



pop <- population("pop1", time = 1, N = N, map = map, competition = 10)
model <- compile_model(pop, generation_time = 1, simulation_length = 1000,
                       resolution = 0.1, mating = 4 , dispersal = 4)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])

p4 <- ggplot(res = 1500) + geom_sf(data = map) + geom_sf(data = nodes, col = rainbow(16)[k])
plot(p4)
gs4 <- fviz_gap_stat(gap_stat, linecolor = 'brown1')



pop <- population("pop1", time = 1, N = N, map = map, competition = 2.25)
model <- compile_model(pop, generation_time = 1, simulation_length = 1000,
                       resolution = 0.1, mating = 0.3 , dispersal = 0.3)

ts <- slim(model, sequence_length = 1, recombination_rate = 0, random_seed = 18181010, burnin = 1000, load = T)
nodes <- ts_nodes(ts) %>% filter(time == max(time)) %>% distinct(ind_id, .keep_all = TRUE)
points <- st_coordinates(nodes)
set.seed(1810)
gap_stat <- clusGap(points, FUN = hcut,nstart = 20, K.max = 15, B = 2, hc_method = 'average',
                    SE.factor = 1,)
k <-  maxSE(f         = gap_stat$Tab[,"gap"],
            SE.f      = gap_stat$Tab[,"SE.sim"])
km <- hcut(points, k, hc_func = 'hclust', hc_method = 'average')

p5 <- ggplot(res = 1500) + geom_sf(data = map) + geom_sf(data = nodes, col = ifelse(km$cl == 1, 'aquamarine3', ifelse(km$cl == 2, 'coral', 'orchid2')))
plot(p5)
gs5 <- fviz_gap_stat(gap_stat, linecolor = 'orchid')
k

ggpubr::ggarrange(p1, p3, p4, p5, gs1,gs3, gs4, gs5, ncol = 4, nrow = 2)

ggpubr::ggarrange(p4, p1, p5, p3, ncol = 4, nrow = 1)
ggpubr::ggarrange(gs4,gs1, gs5, gs3, ncol = 2, nrow = 2)

