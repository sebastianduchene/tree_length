library(cluster)

library(ape)
library(doParallel)

dat_raw <- as.matrix(read.table('topo_dist_mat.txt', head = T, as.is = T))
dat_raw[upper.tri(dat_raw, diag = T)] <- 0
dat_raw[1:ncol(dat_raw), 1:nrow(dat_raw)] <- as.numeric(dat_raw)

dat_dist <- as.dist(dat_raw)


mds_dat <- cmdscale(dat_dist, k = 2)

rm(list = c('dat_raw', 'dat_dist'))


write.table(mds_dat, file = 'mds_topo.txt', row.names = T)

### functions

clus_fun <- function(x, kmax) sapply(2:(kmax), function(i) clara(x, k = i)$silinfo$avg.width )
clus_fun_par <- function(x, kmax) foreach(i = 2:(kmax), .combine = cbind) %dopar% cluster::clara(x, k = i)$silinfo$avg.width


get_boot_rep <- function(x){
  boot_mat <- cbind(runif(nrow(x), min(x[, 1]), max(x[, 1])), runif(nrow(x), min(x[, 2]), max(x[, 2])))
}

###
if(T){
cl <- makeCluster(10)
registerDoParallel(cl)
true_dat_clus <- clus_fun_par(mds_dat, 50)
stopCluster(cl)


boot_dat_clus <- list()

cl <- makeCluster(10)
registerDoParallel(cl)
for(i in 1:50){
      print(paste('boot_rep' , i))
      boot_dat_temp <- get_boot_rep(mds_dat)
      boot_dat_clus[[i]] <- clus_fun_par(boot_dat_temp, 50)
}


stopCluster(cl)

}

pdf('topo_clusters.pdf')
par(mfrow = c(2, 1))
plot(2:50, as.numeric(true_dat_clus), type = 'l', lwd = 2, col = 'red', ylim = c(0.3, 0.6), ylab = 'Average silhouette width', xlab = 'Number of topology clusters')
for(i in 1:length(boot_dat_clus)){
      points(jitter(2:50), as.numeric(boot_dat_clus[[i]]), pch = 20, col = 'blue')

}

clu1 <- pam(mds_dat, k = 10)

plot(mds_dat, pch = 20, col = rainbow(10)[clu1$clustering], xlab = 'MDS 1', ylab = 'MDS 2')

dev.off()



write.table(t(t(clu1$clustering)), file = 'clust_k_10.txt', row.names = T)
write.table(clu1$clusinfo, file = 'clust_info_k_10.txt', row.names = F)

