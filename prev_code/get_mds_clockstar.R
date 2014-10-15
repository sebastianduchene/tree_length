library(cluster)
library(ape)
library(doParallel)

dat_raw <- as.matrix(read.table('sbsd.txt', head = T, as.is = T))
dat_raw[upper.tri(dat_raw, diag = T)] <- 0
dat_raw[1:ncol(dat_raw), 1:nrow(dat_raw)] <- as.numeric(dat_raw)

dat_dist <- as.dist(dat_raw)
mds_dat <- cmdscale(dat_dist, k = 3)




### functions

clus_fun <- function(x, kmax) sapply(2:(kmax), function(i) clara(x, k = i)$silinfo$avg.width )
clus_fun_par <- function(x, kmax) foreach(i = 2:(kmax), .combine = cbind) %dopar% cluster::clara(x, k = i)$silinfo$avg.width

get_boot_rep <- function(x){
  boot_mat <- cbind(runif(nrow(x), min(x[, 1]), max(x[, 1])), runif(nrow(x), min(x[, 2]), max(x[, 2])), runif(nrow(x), min(x[, 3]), max(x[, 3])))
}

###

cl <- makeCluster(4)
registerDoParallel(cl)
true_dat_clus <- clus_fun_par(mds_dat, 10)
stopCluster(cl)

if(T){
boot_dat_clus <- list()

cl <- makeCluster(4)
registerDoParallel(cl)
for(i in 1:50){
      print(paste('boot_rep' , i))
      boot_dat_temp <- get_boot_rep(mds_dat)
      boot_dat_clus[[i]] <- clus_fun_par(boot_dat_temp, 10)
}


stopCluster(cl)
}

#stop('clustering up to here')

par(mfrow = c(2, 1))
plot(2:10, as.numeric(true_dat_clus), type = 'l', lwd = 2, col = 'red', ylim = c(0, 0.7), ylab = 'Average silhouette width', xlab = 'Number of pacemakers')
for(i in 1:length(boot_dat_clus)){
      points(jitter(2:10), as.numeric(boot_dat_clus[[i]]), pch = 20, col = rgb(0, 0, 0.5, 0.2))
#      points(jitter(2:10), as.numeric(boot_dat_clus[[i]]), pch = 20, col = 1)

}



clu1 <- pam(mds_dat, k = 4)

plot(mds_dat, pch = 20, col = clu1$clustering, xlab = 'MDS 1', ylab = 'MDS 2')

write.table(clu1$clustering, file = 'clustering_k_3.txt')
write.table(clu1$clusinfo, file = 'clusinfo_k_3.txt')