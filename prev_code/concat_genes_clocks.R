library(ape)


clock_table <- read.table(grep('clustering_k_.+[.]txt', dir(), value = T))


concat_genes <- list()
for(i in 1:nrow(unique(clock_table))){
      genes_names <- rownames(clock_table)[clock_table == i]
      concat_genes[[i]] <- read.dna(genes_names[1], format = 'fasta')
      for(k in 2:length(genes_names)){
      	    concat_genes[[i]] <- cbind(concat_genes[[i]], read.dna(genes_names[k], format = 'fasta'))
      }
      write.dna(concat_genes[[i]], file = 'concat_clocks.phy', append = T, nbcol = -1, colsep = '')
}


raw_lines <- readLines('concat_clocks.phy')
raw_lines[2:length(raw_lines)] <- gsub(' {1,}', '          ', raw_lines[2:length(raw_lines)])
cat(raw_lines, sep = '\n', file = 'concat_clocks.phy')

tree_concat <- read.tree(grep('nj.+concat.+tree', dir(), value = T))

tree_concat$edge.length <- NULL
tree_concat <- root(tree_concat, outgroup = 's1', resolve.root = T)
cat(nrow(concat_genes[[1]]), file = 'concat_tree.phy', sep = '\n')
write.tree(tree_concat, file = 'concat_tree.phy', append = T)