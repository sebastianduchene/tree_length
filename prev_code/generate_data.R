# Data simulation for genome dating

## Tree topology 

### use 

library(phangorn)

trees_topo <- read.tree('simtree1.tre')

for(i in 1:length(trees_topo)){
  trees_topo[[i]] <- chronopl(trees_topo[[i]], lambda = 5, age.min = 100)
  print(max(branching.times(trees_topo[[i]])))

}

write.tree(trees_topo, file = 'sim_chrono.tre')

set.seed(123456)
pm1 <- rlnorm(38, -6.9, 0.2)

set.seed(654421)
pm2 <- rlnorm(38, -6.9, 0.2)

set.seed(222222)
pm3 <- rlnorm(38, -6.9, 0.2)

set.seed(33333)
pm4 <- rlnorm(38, -6.9, 0.2)

set.seed(666666)
pm5 <- rlnorm(38, -6.9, 0.2)


# important q's: 
  	    	 #can we recover the clusters?
  	    	 # How far are we from the gene tree
var_sites_mat <- matrix(NA, 510, 2)
var_temp <- 0

#tree 1. Large wide cluster with slow pacemakers
      tree_temp <- trees_topo[[1]]
      tree_temp$edge.length <- tree_temp$edge.length * pm1 / 2
      for(k in 1:50){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_', 1, '_pm_1_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')
	    
	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_', 1, '_pm_1_', k, '.fasta'), length(seg.sites(sim_temp_1)))

      }

      tree_temp <- trees_topo[[1]]
      tree_temp$edge.length <- tree_temp$edge.length * pm2 / 2
      for(k in 1:40){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_', 1, '_pm_2_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')

	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_', 1, '_pm_2_', k, '.fasta'), length(seg.sites(sim_temp_1)))

      }

      tree_temp <- trees_topo[[1]]
      tree_temp$edge.length <- tree_temp$edge.length * pm3 / 1.5
      for(k in 1:30){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_', 1, '_pm_3_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')

	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_', 1, '_pm_3_', k, '.fasta'), length(seg.sites(sim_temp_1)))

	}

#tree 2. Large narrow cluster with fast pacemakers
      tree_temp <- trees_topo[[2]]
      tree_temp$edge.length <- tree_temp$edge.length * pm1  * 1.5
      for(k in 1:40){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_', 2, '_pm_1_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')

	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_', 2, '_pm_1_', k, '.fasta'), length(seg.sites(sim_temp_1)))

      }

      tree_temp <- trees_topo[[2]]
      tree_temp$edge.length <- tree_temp$edge.length * pm2  * 1.5
      for(k in 1:30){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_', 2, '_pm_2_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')

	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_', 2, '_pm_2_', k, '.fasta'), length(seg.sites(sim_temp_1)))

      }

      tree_temp <- trees_topo[[2]]
      tree_temp$edge.length <- tree_temp$edge.length * pm3 * 1.5
      for(k in 1:30){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_', 2, '_pm_3_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')

	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_', 2, '_pm_3_', k, '.fasta'), length(seg.sites(sim_temp_1)))

	}


if(T){
#tree 3. Large narrow cluster with fast pacemakers
      tree_temp <- trees_topo[[3]]
      tree_temp$edge.length <- tree_temp$edge.length * pm1  * 1.5 
      for(k in 1:40){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_', 3, '_pm_1_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')

	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_', 3, '_pm_1_', k, '.fasta'), length(seg.sites(sim_temp_1)))

      }

      tree_temp <- trees_topo[[3]]
      tree_temp$edge.length <- tree_temp$edge.length * pm2  * 1.5
      for(k in 1:30){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_', 3, '_pm_2_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')

	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_', 3, '_pm_2_', k, '.fasta'), length(seg.sites(sim_temp_1)))

      }

      tree_temp <- trees_topo[[3]]
      tree_temp$edge.length <- tree_temp$edge.length * pm3 
      for(k in 1:30){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_', 3, '_pm_3_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')

	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_', 3, '_pm_3_', k, '.fasta'), length(seg.sites(sim_temp_1)))

	}
}

if(T){
#tree 4. Large narrow cluster with fast pacemakers
      tree_temp <- trees_topo[[4]]
      tree_temp$edge.length <- tree_temp$edge.length * pm1  / 2
      for(k in 1:40){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_', 4, '_pm_1_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')

	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_', 4, '_pm_1_', k, '.fasta'), length(seg.sites(sim_temp_1)))

      }

      tree_temp <- trees_topo[[4]]
      tree_temp$edge.length <- tree_temp$edge.length * pm2  / 2
      for(k in 1:30){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_', 4, '_pm_2_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')

	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_', 4, '_pm_2_', k, '.fasta'), length(seg.sites(sim_temp_1)))

      }

      tree_temp <- trees_topo[[4]]
      tree_temp$edge.length <- tree_temp$edge.length * pm3 / 2
      for(k in 1:30){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_', 4, '_pm_3_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')

	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_', 4, '_pm_3_', k, '.fasta'), length(seg.sites(sim_temp_1)))

	}
}



#tree 5. Large narrow cluster with fast pacemakers
      tree_temp <- trees_topo[[5]]
      tree_temp$edge.length <- tree_temp$edge.length * pm1  * 2.5
      for(k in 1:40){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_', 5, '_pm_1_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')

	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_', 5, '_pm_1_', k, '.fasta'), length(seg.sites(sim_temp_1)))

      }

      tree_temp <- trees_topo[[5]]
      tree_temp$edge.length <- tree_temp$edge.length * pm2 * 2.5
      for(k in 1:30){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_', 5, '_pm_2_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')

	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_', 5, '_pm_2_', k, '.fasta'), length(seg.sites(sim_temp_1)))

      }

      tree_temp <- trees_topo[[5]]
      tree_temp$edge.length <- tree_temp$edge.length * pm3 * 3.5
      for(k in 1:20){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_', 5, '_pm_3_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')

	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_', 5, '_pm_3_', k, '.fasta'), length(seg.sites(sim_temp_1)))


	}


write.table(var_sites_mat, file = 'var_stes.txt', row.names = F)




 
=======

var_sites_mat <- matrix(NA, 500, 2)
var_temp <- 0

#tree 1. Large wide cluster with low rates
      tree_temp <- trees_topo[[1]]
      tree_temp$edge.length <- tree_temp$edge.length * pm1 / 3
      for(k in 1:150){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_1_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')
	    
	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_1_', k, '.fasta'), length(seg.sites(sim_temp_1)))

      }

#tree 2. Large wide cluster with low rates
      tree_temp <- trees_topo[[2]]
      tree_temp$edge.length <- tree_temp$edge.length * pm1 / 3
      for(k in 1:100){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_2_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')
	    
	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_2_', k, '.fasta'), length(seg.sites(sim_temp_1)))

      }

#tree 3. Large wide cluster with low rates
      tree_temp <- trees_topo[[3]]
      tree_temp$edge.length <- tree_temp$edge.length * pm1 / 3
      for(k in 1:100){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_3_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')
	    
	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_3_', k, '.fasta'), length(seg.sites(sim_temp_1)))

      }

#tree 4. Large wide cluster with low rates
      tree_temp <- trees_topo[[4]]
      tree_temp$edge.length <- tree_temp$edge.length * pm1 / 3
      for(k in 1:100){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_4_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')
	    
	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_4_', k, '.fasta'), length(seg.sites(sim_temp_1)))

      }

#tree 5. Large wide cluster with low rates
      tree_temp <- trees_topo[[5]]
      tree_temp$edge.length <- tree_temp$edge.length * pm1
      for(k in 1:50){
      	    sim_temp_1 <- as.DNAbin(simSeq(tree_temp, l = 1000))
	    write.dna(sim_temp_1, file = paste0('tr_5_', k, '.fasta'), format = 'fasta', nbcol = -1, colsep = '')
	    
	    var_temp <- var_temp + 1
	    var_sites_mat[var_temp, ] <- c(paste0('tr_5_', k, '.fasta'), length(seg.sites(sim_temp_1)))

      }


write.table(var_sites_mat, file = 'var_stes.txt', row.names = F)

>>>>>>> 3267e3604a56094491099b7f02411c212f75a6de
