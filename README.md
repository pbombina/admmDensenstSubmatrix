```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "vignettes/icon-",
  out.width = "100%"
)
```


# admmDensenstSubmatrix

# Introduction
This is the `R`-package accompanying the paper ([Convex optimization for the densest subgraph and densest submatrix problems](https://github.com/bpames/Densest-Submatrix-Paper/blob/master/Manuscript/dsm-arxiv2019.pdf).

The problem of identifying a dense submatrix is a fundamental problem in the  analysis of matrix structure and complex networks. This package provides tools for identifying the densest submatrix of a given graph using first-order optimization methods.

See the tutorials below to get started.

# The densest submatrix problem
Let ![](https://latex.codecogs.com/gif.latex?%5BM%5D%20%3D%20%5C%7B1%2C2%2C%5Cdots%2C%20M%5C%7D)for each positive integer ![](https://latex.codecogs.com/gif.latex?M).
Given a matrix ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BA%7D%20%5Cin%20R%5E%7BM%5Ctimes%20N%7D), the densest ![](https://latex.codecogs.com/gif.latex?%24m%5Ctimes%20n%24)-submatrix problem seeks subsets ![](https://latex.codecogs.com/gif.latex?%5Cbar%20U%20%5Csubseteq%20%7B%5BM%5D%7D) and ![](https://latex.codecogs.com/gif.latex?%5Cbar%20V%20%5Csubseteq%20%7B%5BN%5D%7D) of cardinality ![](
https://latex.codecogs.com/gif.latex?%7C%5Cbar%20U%7C%3Dm) and ![](https://latex.codecogs.com/gif.latex?%7C%5Cbar%20V%7C%20%3D%20n), respectively,
such that the submatrix ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BA%7D%7B%5B%5Cbar%20U%2C%20%5Cbar%20V%5D%7D) with rows index by ![](https://latex.codecogs.com/gif.latex?%5Cbar%20U) and columns indexed by ![](https://latex.codecogs.com/gif.latex?%5Cbar%20V)
contains the maximum number of nonzero entries. That is, the densest ![](https://latex.codecogs.com/gif.latex?%24m%5Ctimes%20n%24)-submatrix problem seeks the densest
![](https://latex.codecogs.com/gif.latex?%24m%5Ctimes%20n%24)-submatrix of ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BA%7D).

The densest ![](https://latex.codecogs.com/gif.latex?m%5Ctimes%20n)-submatrix problem can be formulated as:

![](https://latex.codecogs.com/gif.latex?%5Cmin_%7B%5Cmathbf%7BX%7D%2C%20%5Cmathbf%7BY%7D%20%5Cin%20%7B%20%5C%7B0%2C1%5C%7D%7D%5E%7BM%5Ctimes%20N%7D%20%7D%20%7B%5Ctr%28%5Cmathbf%7BY%7D%20%5Cmathbf%7Be%7D%20%5Cmathbf%7Be%7D%5ET%29%3A%20%5Cmathrm%7BP%7D_%7B%5COmega%7D%28%5Cmathbf%7BX%7D-%5Cmathbf%7BY%7D%29%20%3D%20%5Cmathbf%7B0%7D%2C%20tr%28%5Cmathbf%7BX%7D%20%5Cmathbf%7Be%7D%20%5Cmathbf%7Be%7D%5ET%29%20%3D%20mn%2C%20rank%20%28%5Cmathbf%7BX%7D%29%20%3D%201%20%7D)
  
where

* ![](https://latex.codecogs.com/gif.latex?%5Cmathrm%7BP%7D_%7B%5COmega%7D)is the projection onto the index set of zero entries of matrix ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%20A);

* ![](https://latex.codecogs.com/gif.latex?tr) is the matrix trace function;

* ![](https://latex.codecogs.com/gif.latex?%5COmega) is the index set of zero entries of matrix ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%20A);

* ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BX%7D) is rank-one matrix with ![](https://latex.codecogs.com/gif.latex?mn) nonzero entries;

* ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BY%7D) is used to count the number of disagreements between ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BA%7D) and ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BX%7D);

* ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7Be%7D) - all-ones vector.

Unfortunately, optimization problems involving rank and binary constraints are intractable in general.

Relaxing the rank constraint with a nuclear norm penalty term,
![](https://latex.codecogs.com/gif.latex?%5C%7C%5Cmathbf%7BX%7D%20%5C%7C_*%20%3D%20%5Csum_%7Bi%3D1%7D%5EN%20%5Csigma_i%28%5Cmathbf%7BX%7D%29), i.e., the sum of the singular values of matrix, and
the binary constraints with box constraints yields the convex problem:

![](https://latex.codecogs.com/gif.latex?%5Cbegin%7Balign*%7D%20%5Cmin%20%5C%3B%20%26%20%5C%7C%5Cmathbf%7BX%7D%20%5C%7C_*%20&plus;%20%5Cgamma%20tr%28%5Cmathbf%7BY%7D%20%5Cmathbf%7Be%7D%20%5Cmathbf%7Be%7D%5ET%29%20%5C%5C%20s.t.%20%5C%3B%20%26%20tr%28%5Cmathbf%7BX%7D%20%5Cmathbf%7Be%7D%20%5Cmathbf%7Be%7D%5ET%29%20%3D%20mn%2C%20%5C%5C%20%26%20%5Cmathrm%7BP%7D_%5COmega%28%5Cmathbf%7BX%7D%20-%20%5Cmathbf%7BY%7D%29%20%3D%20%5Cmathbf%7B0%7D%2C%20%5C%5C%20%26%20%5Cmathbf%7BY%7D%20%5Cge%20%5Cmathbf%7B0%7D%2C%20%5C%5C%20%26%20%5Cmathbf%7B0%7D%20%5Cle%20%5Cmathbf%7BX%7D%20%5Cle%20%5Cmathbf%7Be%7D%20%5Cmathbf%7Be%7D%5ET%2C%20%5Cend%7Balign*%7D)

where ![](https://latex.codecogs.com/gif.latex?%5Cgamma%20%3E0) is a regularization parameter chosen to tune between the two objectives.

It can be shown that the relaxation is exact when binary matrix ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BA%7D) contains a single, relatively large dense ![](https://latex.codecogs.com/gif.latex?m%5Ctimes%20n) block. For more information, see ([Convex optimization for the densest subgraph and densest submatrix problems](https://github.com/bpames/Densest-Submatrix-Paper/blob/master/Manuscript/dsm-arxiv2019.pdf))

# Alternating Direction Method of multipliers for densest submatrix problem
The alternating direction method of multipliers (ADMM) has been succesfully used in a broad spectrum of applications. The ADMM solves convex optimization problems with composite objective functions subject to equality constraints.  

We direct the reader to Prof. Stephen Boyd’s website
([ADMM](http://stanford.edu/~boyd/papers/admm_distr_stats.html)) for a more thorough discussion of the ADMM.

To apply ADMM to our problem, we introduce artificial variables ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BQ%7D), ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BW%7D) and ![](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BZ%7D) to obtain the equivalent optimization problem:

![](https://latex.codecogs.com/gif.latex?%5Cbegin%7Balign*%7D%20%5Cmin%20%5C%3B%20%26%20%5C%7C%5Cmathbf%7BX%7D%20%5C%7C_*%20&plus;%20%5Cgamma%20%5C%7C%5Cmathbf%7BY%7D%20%5C%7C_1%20&plus;%7B1%7D_%7B%5COmega_Q%7D%28%5Cmathbf%7BQ%7D%29&plus;%7B1%7D_%7B%5COmega_W%7D%28%5Cmathbf%7BW%7D%29&plus;%7B1%7D_%7B%5COmega_Z%7D%28%5Cmathbf%7BZ%7D%29%5C%5C%20s.t.%20%5C%3B%20%26%20%5Cmathbf%7BX%7D-%5Cmathbf%7BY%7D%3D%5Cmathbf%7BQ%7D%2C%5Cmathbf%7BX%7D%3D%5Cmathbf%7BW%7D%2C%20%5Cmathbf%7BX%7D%3D%5Cmathbf%7BZ%7D%20%5Cend%7Balign*%7D)

where

* ![](https://latex.codecogs.com/gif.latex?%5COmega_Q%20%3D%20%5C%7B%5C%2C%20%5Cmathbf%7BQ%7D%5Cin%20R%5E%7BM%5Ctimes%20N%7D%20%5Cmid%20P_%7B%5Ctilde%7BN%7D%7D%28Q%29%3D0%20%5C%2C%20%5C%7D),
* ![](https://latex.codecogs.com/gif.latex?%5COmega_W%20%3D%5C%7B%5C%2C%20%5Cmathbf%7BW%7D%5Cin%20R%5E%7BM%5Ctimes%20N%7D%20%5Cmid%20%5Cmathbf%7Be%7D%5ET%5Cmathbf%7BW%7D%20%5Cmathbf%7Be%7D%3Dmn%20%5C%2C%20%5C%7D),
* ![](https://latex.codecogs.com/gif.latex?%5COmega_Z%20%3D%5C%7B%5C%2C%20%5Cmathbf%7BZ%7D%5Cin%20R%5E%7BM%5Ctimes%20N%7D%20%5Cmid%20%7BZ%7D_%7Bij%7D%5Cleq%201%20%5Cforall%20%28i%2Cj%29%5Cin%20M%5Ctimes%20N%20%5C%2C%20%5C%7D).

Here ![](https://latex.codecogs.com/gif.latex?%7B1%7D_%7BS%7D%3A%20R%5E%7BM%5Ctimes%20M%7D%20%5Crightarrow%20%5Cleft%20%5C%7B0%2C&plus;%5Cinfty%20%5Cright%20%5C%7D)  is the indicator function of the set ![](https://latex.codecogs.com/gif.latex?S%20%5Csubseteq%20R%5E%7BM%5Ctimes%20N%7D),
such that
![](https://latex.codecogs.com/gif.latex?%7B1%7D_S%28%5Cmathbf%7BX%7D%29%3D0%24%20if%20%24%5Cmathbf%7BX%7D%5Cin%20S), and ![](https://latex.codecogs.com/gif.latex?&plus;%5Cinfty) otherwise.

Since our objective function is separable, we iteratively solve this optimization program using the ADMM.
The basic idea is to rotate through 3 steps:

1. minimize the augmented Lagrangian over primal variables,
2. update dual variables usng the updated primal variables,
3. calculate primal and dual residuals.

Interested readers are referred to ([Convex optimization for the densest subgraph and densest submatrix problems](https://github.com/bpames/Densest-Submatrix-Paper/blob/master/Manuscript/dsm-arxiv2019.pdf)). We include a summary of the algorithm below.

![](https://github.com/pbombina/admmDensenstSubmatrix/blob/master/vignettes/ALG.png?raw=true)

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

