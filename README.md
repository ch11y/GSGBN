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
1.  The package does not contain the process cross validation.
2.  For the details on how to run the simulated and real datasets, please refer to the README files in the directory of the data.
