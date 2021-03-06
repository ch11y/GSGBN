The tool kit contains GSGBN's software for ICDM14 paper *"Learning Sparse Gaussian Bayesian Networks by Variable Grouping"*.

#INSTALL
Exract the package, then use make to compile the source code.
-  $ git clone https://github.com/ch11y/GSGBN/
-  $ cd GSGBN
-  $ make

Find the excutive file in ./bin/.
#RUN
##Usage:
    GSGBN -i input.txt -l1 1.0 -l2 0.1 -num 10 -thr 0.1
##Options:
    -i,   an m * n matrix with m samples and n variables;
    -l1,  the regularization parameter lambda1 for sparsity;
    -l2,  the regularization parameter lambda2*n is used for grouping;
    -num, number of rounds for enumerating DAGs;
    -thr, threshold for filtering
##*Notes*
***The package does not contain the process of cross validation.***

##Contact:
Jie Yang (Email: jyang2@cs.hku.hk), Department of Computer Science, The University of Hong Kong.
#Reference
Jie Yang, Henry C.M. Leung, S.M. Yiu, Yunpeng Cai, Francis Y.L. Chin, *Learning Sparse Gaussian Bayesian Networks by Variable Grouping*, The IEEE International Conference on Data Mining (ICDM 2014).
