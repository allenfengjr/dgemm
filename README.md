# DGemm Report

Hao Feng haofeng@iu.edu

## Problem Restate


Calculating matrix-multiplication on one computer is a challenging work when matrices are very big. Sometimes one computer is even not able to store the whole matrix. In that case, I implement Cannon's Algorithm to calculate matrix-matrix production.



## Final Parallel Design

My final parallel design includes two level parallelization. In Cannon's Algorithm, there are many block



The first one is parallel using MPI. It can both run on inter-nodes and intra-node. Using MPI can enhance the scalability of the entire system. Especially on weak scaling.



The second one is parallel using OpenMP. It can only run on intra-node because it use shared memory. Using OpenMP can enhance strong scaling on each node.



## Platform Information

I use BigRed 200 as my experiment platform.



**Spec of Platform**

Big Red 200 features 640 compute nodes, each equipped with 256 GB of memory and two 64-core, 2.25 GHz, 225-watt AMD EPYC 7742 processors.



## Code

My code has been uploaded to [DGemm](https://github.com/allenfengjr/dgemm).
