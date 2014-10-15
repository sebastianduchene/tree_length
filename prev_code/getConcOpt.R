require(phangorn)

getConcOpt <- function(fastastring, out_files = 'topo_clus_1'){

	   genes <- list(read.dna(fastastring[1], format = "fasta"))
	   concdat <- genes[[1]]
	   for(i in 2:length(fastastring)){
	   	 genes[[i]] <- read.dna(fastastring[i], format = "fasta")
	   	 concdat <- cbind(concdat, genes[[i]])
	   }
	   
	   conctr <- optim.pml(pml(nj(dist.dna(concdat)), as.phyDat(concdat)))$tree
	   gentrs <- list()
	   for(i in 1:length(genes)){
		 gentrs[[i]] <- optim.pml(pml(conctr, as.phyDat(genes[[i]])))$tree
		 
	   }

	   class(gentrs) <- "multiPhylo"
	   names(gentrs) <- gsub('^.+/', '', fastastring)
	   print(gentrs)
	   
	   write.tree(conctr, file = paste0(out_files, '_concat.tree'))
	   write.tree(gentrs, file = paste0(out_files, '_gts.trees'), tree.names = T)


}

fasta_files_1 = grep('tr_1_*', dir(pattern = 'fasta', rec = T), value = T)
fasta_files_5 = grep('tr_5_*', dir(pattern = 'fasta', rec = T), value = T)

getConcOpt(fasta_files_1, out_files = 'topo_clus_1')

getConcOpt(fasta_files_5, out_files = 'topo_clus_5')

