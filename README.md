Neighbor Join Tree based sample selection with Maximal Phylogenetic Diversity
-----------------------------------------

* Obtain

    make
    ./main infile.txt outfile.txt result
    
* Use

   ./main [file1] [file2] [output_file]

* Contact

    This program is jointly developed by Xiaowei Zhan (zhanxw@gmail.com) and Peng Zhang (behappyzp@gmail.com). Comments and feedbacks are welcomed.

This program is associated with publication *Genotype Imputation Reference Panel Selection Using Maximal Phylogenetic Diversity* by Peng Zhang et al. (2013) Genetics. 195(2): 319-330 [link](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3781962/). Note because this hasn't been developed for general use perpuse so the users may need to adapt the code for their own purpose and also need to follow the exact format as it showed in the example files.

The code shown here is for the implementation for maximal phylogenetic diversity algorithm in C++. The program takes two input files, one showed in the format of *infile.txt*, which is a pair-wise Hamming distance matrix between individuals. The Hamming distances are computed by counting sequence differences between a pair of individuals across the region of comparison. Note that individuals need to be named as showed here, ind1, ind2, ind3, ..., and so on. The first line is the sample size, and should have spaces before the number as showed in the file.

The second file should be in the format of *outfile.txt*, each line indicates the distance between two nodes on the neighbor-joining (NJ) tree. The NJ tree is constructed using the implementation by PHYLIP package (http://evolution.genetics.washington.edu/phylip/doc/neighbor.html), where it takes the input file with format shown in *infile.txt* and generates a NJ tree. The *outfile.txt* is part of program output file. 

[![Build Status](https://travis-ci.org/zhanxw/SPD.svg)](https://travis-ci.org/zhanxw/SPD)
