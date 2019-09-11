# admmDensenstSubmatrix  <img src="vignettes/icon.png" align="right" height=170/>


# Introduction
This is the `R`-package accompanying the paper ([Convex optimization for the densest subgraph and densest submatrix problems](https://github.com/bpames/Densest-Submatrix-Paper/blob/master/Manuscript/dsm-arxiv2019.pdf).

The problem of identifying a dense submatrix is a fundamental problem in the  analysis of matrix structure and complex networks. This package provides tools for identifying the densest submatrix of a given graph using first-order optimization methods.

See the tutorials below to get started.

## Installation

```r
#Install the development version from GitHub:
# install.packages("remotes")
remotes::install_github("pbombina/admmDensenstSubmatrix")

```
To also build the vignettes use:

```r
#install.packages("remotes")
remotes::install_github("pbombina/admmDensenstSubmatrix", dependencies = TRUE,
                         build_vignettes = TRUE)

```

## Overview
In this study we consider the densest submtarix problem which seeks to find the densest submatrix of the given size.

We seek the submatrix of the desired size with maximum number of nonzero elements. We proposed a new convex relaxation to solve this problem which based on nuclear norm relaxation. Our relaxation correctly identifies the densestsubmatrix of the fixed size in random matrices if the entries within this submatrix are significantly more likely to be nonzero than arbitrary entry of the matrix.  

See the paper for more details, in particular regarding the derivation of Alternating Direction Method of Multipliers (ADMM) for densest submatrix problem.


## Contents of this repository
`R` scripts are organized into the following subdirectories. 

- [data_preparation](data_preparation/): preparation of benchmark data files
- [ensemble_clustering](ensemble_clustering/): run and evaluate ensemble clustering
- [evaluate_results](evaluate_results/): scripts to evaluate results from all methods
- [helpers](helpers/): helper functions
- [plots_and_tables](plots_and_tables/): generate plots and tables of results
- [range_k](range_k/): run and evaluate FlowSOM over range of values k (number of clusters)
- [run_methods](run_methods/): scripts to run all methods (or instructions to run graphical interfaces, where required)
- [stability_analysis](stability_analysis/): run and evaluate methods for stability analysis

Supplementary files from the published paper are included in the following directory:

- [supplementary_files](supplementary_files/): supplementary files from paper (latest version: November 18, 2016)

R scripts and summary reports for updated results are included in the following directory:

- [updates](updates/): updated results for new clustering methods or new reference data sets


# Examples
We test this package on two different types of data: first, using random matrices sampled from the planted dense $m \times n$ submtarix model and, second, real-world collaboration and communication networks.

## Random matrices
We first generate a random matrix with noise obscuring the planted submatrix using the function ``plantedsubmatrix``. and then call the function ``densub`` to recover the planted submatrix.

```R
# Initialize problem size and densities
# You can play around with these parameters
M <- 100 #number of rows of sampled matrix
N <- 200 #number of columns of sampled matrix
m <- 50 #number of rows of dense submatrix
n <- 40 #number of columns of dense submatrix
p <- 0.25 # noise density
q <- 0.85 #in-group density

#Make binary matrix with mn-submatrix
random<-plantedsubmatrix(M = M, N = N,m = m,n = n,p = p,q = q)
```

After generating the structure `random` containing the random matrix with desired planted structure, we can visually represent the matrix and planted submatrix as two-tone images, where dark pixels correspond to nonzero entries, and light pixels correspond to zero entries, using the following commands.

```R

# Plot sampled G and matrix representations.
image(random$sampled_matrix, useRaster = TRUE, axes = FALSE, main = "Matrix A")
image(random$dense_submatrix, useRaster = TRUE, axes = FALSE, main = "Matrix X0")
image(random$disagreements, useRaster = TRUE, axes = FALSE, main = "Matrix Y0")
```

Tne vizualization of the randomly generated matrix ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BA%7D) helps us to understand its structure. It is clear that ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BA%7D) contains a dense ![](https://latex.codecogs.com/gif.latex?50%20%5Ctimes%2040) block (in the bottom left corner).

![Visual representation of randomly generated ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BA%7D)](https://github.com/pbombina/admmDensenstSubmatrix/blob/master/vignettes/Rplot.jpeg?raw=true)

We can remove all noise and isolate an image of a rank-one matrix ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BX0%7D) with ![](https://latex.codecogs.com/gif.latex?mn) nonzero entries.

![Visual representation of dense submatrix](https://github.com/pbombina/admmDensenstSubmatrix/blob/master/vignettes/Rplot01.jpeg?raw=true)


Then we vizualize matrix ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BY0%7D) to see the number of disagreements between original matrix ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BA%7D) and ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BX0%7D).

![Disagreement between $\mathbf{A}$ and $\mathbf{X_0}$](https://github.com/pbombina/admmDensenstSubmatrix/blob/master/vignettes/Rplot02.jpeg?raw=true)

We call the ADMM solver and visualize the output using the following commands.


```R
#Call ADMM solver
admm <- densub(G = random$sampled_matrix, m = m, n = n, tau = 0.35, gamma = 6/(sqrt(m*n)*(q-p)), opt_tol = 1.0e-4,maxiter = 500, quiet = TRUE)


#Plot results
image(admm$X, useRaster = TRUE, axes = FALSE, main = "Matrix X")
image(admm$Y, useRaster = TRUE, axes = FALSE, main = "Matrix Y")


```


The ADMM solver returns the optimal solutions ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BX%7D) and ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BY%7D). It must be noted that matrices ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BX%7D) and ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BY%7D) are identical to the actual structures of ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BX0%7D) and ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BY0%7D). The planted submatrix is recovered.

![Optimal solution \mathbf{X}](https://github.com/pbombina/admmDensenstSubmatrix/blob/master/vignettes/Rplot03.jpeg?raw=true)


![Optimal Solution \mathbf{Y}](https://github.com/pbombina/admmDensenstSubmatrix/blob/master/vignettes/Rplot04.jpeg?raw=true)


## Collaboration Network
The following is a simple example on how one could use the package to analyze the collaboration network found in the JAZZ dataset. It is known that this network contains a cluster of 100 musicians which performed together.

![JAZZ Network](https://github.com/pbombina/admmDensenstSubmatrix/blob/master/vignettes/0001.jpg?raw=true)

We have already prepared dataset to work with. More details can be found in the provided file `JAZZ_IN_R.R.`

```R
#Load dataset
load(file = "JAZZ.RData")

#Initialize problem size and densities
G <- new #define matrix G equivalent to JAZZ dataset 
m <- 100 #clique size or the number of rows of the dense submatrix 
n <- 100 #clique size of the number of columns of the dense sumbatrix
ta <- 0.85 #regularization parameter
opt_tol <- 1.0e-2 #optimal tolerance
verbose <- 1
maxiter <- 2000 #number of iterations
gamma <- 8/n #regularization parameter



#call ADMM solver
admm <- densub(G = G, m = m, n = n, tau = tau, gamma = gamma, opt_tol = opt_tol, maxiter=maxiter, quiet = TRUE) 
# Planted solution X0.
X0 <- matrix(0L, nrow = 198, ncol = 198) #construct rank-one matrix X0
X0[1:100,1:100] <- matrix(1L, nrow = 100, ncol = 100)#define dense block

# Planted solution Y0.
Y0 <- matrix(0L, nrow = 198, ncol = 198) #construct matrix for counting disagreements between G and X0
Y0[1:100,1:100] < matrix(1L,nrow = 100,ncol = 1000)-G[1:100,1:100]

#Check primal and dual residuals.
C <- admm$X-X0
a <- norm(C, "F") #Frobenius norm of matrix C 
b <- norm(X0,"F") #Frobenius norm of matrix X0
recovery <- matrix(0L,nrow = 1, ncol = 1)#create recovery matrix

if (a/b^2<opt_tol){
recovery = recovery+1
} else {
  recovery = 0 #Recovery condition 
  }





```

Our algorithm converges to the dense submatrix representing the community of 100 musicians after 50 iterations.     

