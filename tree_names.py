import re
import numpy as np
import pandas as pd

input_tree_file = raw_input('Please drag your tree file here:\n')
#15_taxa_unique.tree
tree_dat = np.loadtxt(re.sub(' ', '', input_tree_file), delimiter = '(', dtype = str)
tree_names = pd.Series(tree_dat[:, 0])
tree_name_counts = tree_names.value_counts()

if np.any(tree_name_counts > 1):
    print "The following trees are duplicated in the file: " %(tree_name_counts.index[tree_name_counts != 1][0])
else: 
    print "There are no duplicate trees in this file"
