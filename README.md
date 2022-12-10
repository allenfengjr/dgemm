# DGemm Report

Hao Feng haofeng@iu.edu

## Problem Restate


Calculating matrix-multiplication on one computer is a challenging work when matrices are very big. Sometimes one computer is even not able to store the whole matrix. In that case, I implement Cannon's Algorithm to calculate matrix-matrix production.



## Final Parallel Design

My final parallel design includes two level parallelization. In Cannon's Algorithm, there are `k*k`blocks, blocks swap sub-matrix with each other, and after one swap, they do a sub-matrix production then update their local result. After `k` times swap, root process gather all result and re-assign the data. For the process of sub-matrix production, we can use parallel algorithm to accelerate it. In this program, I apply parallel-for for this.



The first one is parallelization using MPI. It can both run on inter-nodes and intra-node. Using MPI can enhance the scalability of the entire system. Especially on weak scaling.



The second one is parallelization using OpenMP. It can only run on intra-node because it use shared memory. Using OpenMP can enhance strong scaling on each node.



## Final Parallel Implement

### 1. Program Usage and Limitation

This program needs five arguments.

```
argv[1]: height of a block in A
argv[2]: weight of a block in A
argv[3]: weight of a block in B
argv[4]: the value of alpha
argv[5]: the value of beta
```

The program can accept `A,B,C` as rectangle matrix. It has two limitation. The first one is `#process` must be a square number for spilt. The second is that the width and height must be a multiple of `sqrt(#process)`. If not, all data exchange and matrix production process must consider the problem of data size. I am not able to finish that at this time.

### 2. Algorithm Implement

I would like to present the whole algorithm in three parts, and each part will be accompanied by some code.

**(1) Data Preparation**

In this part, the program spilt the whole matrix into blocks. It needs to mention that we will store those matrices in a buffer, then MPI_Scatter them to avoid a lot of MPI_Send and MPI_Recv.

```c
    int ii = p/k;
    int jj = p%k;
    int A_start = ii*k*(m_size*l_size) + jj*l_size;
    for(int i = 0; i < m_size;++i){
        //copy one row
        std::copy(A.begin()+A_start+i*ll,A.begin()+A_start+i*ll+l_size,A_v.begin()+A_count);
        A_count += l_size;
    }
    MPI_Scatter(A_v.data(),m_size*l_size,MPI_DOUBLE,sub_A.data(),m_size*l_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
```

**(2) Permute Sub-Matrix**

Because I only simply distribute sub-matrices, I still need to initialize them by shifting. I use Descartes Communicator to simplify the process. It is clear for a process to shift a block and map ranks to its logical neighbors.

```c
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,0,&cart_comm);
    MPI_Cart_coords(cart_comm,my_rank,2,coords);
    //init movement, for (i,j) move its sub_A i-times and sub_B j-times
    MPI_Cart_shift(cart_comm, 0, coords[1], &udlr[0], &udlr[1]);
    MPI_Cart_shift(cart_comm, 1, coords[0], &udlr[2], &udlr[3]);
```

**(3) Matrix-Matrix Multiplication**

In this part, the program need to shift `k` times and do `k` times Matrix-Matrix Multiplication. Because I use row-wise matrix, I choose `jki-order` for the loop. I also apply `parallel-for` at the outer loop.

```c
#pragma omp parallel for
    for (int i = 0; i < m_size*n_size; ++i) {
        sub_C[i] *= beta;
    }

    for (int p = 0; p < k; ++p) {
#pragma omp parallel for
        for (int j = 0; j < n_size; ++j) {
            for (int l = 0; l < l_size; ++l) {
                for (int i = 0; i < m_size; ++i) {
                    sub_C[i+j*m_size] += alpha * sub_A[i*l_size+l] * sub_B[l*n_size+j];
                }
            }
        }
        MPI_Barrier(cart_comm);
        MPI_Cart_shift(cart_comm,0,1,&udlr[0],&udlr[1]);
        MPI_Cart_shift(cart_comm,1,1,&udlr[2],&udlr[3]);
        MPI_Sendrecv_replace(sub_A.data(),m_size*l_size,MPI_DOUBLE,udlr[2],3,udlr[3],3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Sendrecv_replace(sub_B.data(),l_size*n_size,MPI_DOUBLE,udlr[0],4,udlr[1],4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
```

After those three process, I gather the result and re-order them, it is like part 1 so I do not explain too much.

## Performance Evaluation

### 1. OpenMP on Intra-node

In this part, I focus on the strong scaling. I use same `#process=64`and `#node=8`. I set all three matrices as square matrix.

| #cpus-per-tasks\Block Size | 2000      | 4000      |
|:--------------------------:|:---------:|:---------:|
| 1                          | 827500ms  | 7836423ms |
| 2                          | 420768ms  | 3925765ms |
| 4                          | 218614ms  | 1993759ms |
| 8                          | 87106.8ms | 862829ms  |
| 16                         | 47572.7ms | 418593ms  |

### 2. MPI on Inter-node

In this part, I focus on both strong scaling and weak scaling. I set `Block size` be 1000, 2000 and `#process` be 4, 9, 16, 25. All of these processes are on different nodes.

| #process\ Block Size | 1000            | 2000           |
|:--------------------:|:---------------:|:--------------:|
| 4                    | 16447.5ms       | **187217ms**   |
| 9                    | 24332.7ms       | ***286500ms*** |
| 16                   | **35075.1ms**   | 384616ms       |
| 25                   | ***44843.9ms*** | 530598ms       |

## 3. Trade-off between OpenMP and MPI

One question I would like to know is that if I have the same number of threads, which is better, using only MPI compared to combining OpenMP and MPI.



In this part, I set `#node`and `#tasks-per-node * #cpus-per-task` fixed, then change `#tasks-per-node` and `#cpus-per-tasks` to find out how to choose a better setup. `#node is 4`

| #tasks-per-node x # cpus-per-tasks\Total Size | 4000      | 8000      |
|:---------------------------------------------:|:---------:|:---------:|
| 1x16                                          | 9731.25ms | 85459.2ms |
| 4x4                                           | 9620.12ms | 97414.1ms |
| 16x1                                          | 9230.48ms | 79832.9ms |



## Conclusion

From the experiment, I can get three conclusions.



The first one is that apply multi-thread parallelization on intra-node has good strong scaling and apply MPI on inter-node has a fair strong scaling(bold pairs in second table). 



The second one is that applying MPI on inter-node has weak scaling. Such as the columns in second table, the time cost is about proportional to the `k`.



The third is that choosing to use OpenMP or MPI on single node is not an easy problem. The data movement cost and potential false sharing issues both have an impact on the choice of the appropriate distribution method.



## Platform Information

I use BigRed 200 as my experiment platform. Big Red 200 features 640 compute nodes, each equipped with 256 GB of memory and two 64-core, 2.25 GHz, 225-watt AMD EPYC 7742 processors.



## Code and Hash

My code has been uploaded to [DGemm](https://github.com/allenfengjr/dgemm).

The Hash is .
