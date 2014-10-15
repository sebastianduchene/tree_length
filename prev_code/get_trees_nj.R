library(ape)


trees_names <- grep('.+fasta', dir(), value = T)


trees_nj_list <- list()
for(i in 1:length(trees_names)){
      dat_temp <- read.dna(trees_names[i], format = 'fasta')
      print(paste('I am estimating tree', trees_names[i]))
      dist_temp <- dist.dna(dat_temp)
      trees_nj_list[[i]] <- nj(dist_temp)
}

      names(trees_nj_list) <- trees_names

      class(trees_nj_list) <- 'multiPhylo'

      write.tree(trees_nj_list, file = 'nj_trees.trees', append = F, tree.names = T)
