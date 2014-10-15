library(ape)

fasta_files <- grep('fasta', dir(), value = T)


clus_150 <- sample(fasta_files, size = 50)

for(i in 1:length(clus_150)){
      system(paste('cp', clus_150[i], 'rand_clus_5_2/'))
}


dat <- read.dna(clus_150[1], format = 'fasta')

for(k in 2:length(clus_150)){
      dat <- cbind(dat, read.dna(clus_150[k], format = 'fasta'))
}


nj_concat_tree <- nj(dist.dna(dat))
write.tree(nj_concat_tree, file = 'rand_clus_5_2/nj_concat_clus_5_2.tree')