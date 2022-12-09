//
// Created by Hao Feng on 11/28/22.
//

#include "iostream"
#include "MMIOReader.hpp"
#include <omp.h>
#include "mpi.h"
#include <chrono>                   // std::chrono::high_resolution_clock::now()
#include <cmath>                    // std::sin, M_PI
#include <numeric>                  // std::transform_reduce
#include "Matrix.hpp"
#include "Vector.hpp"
/*
 * args usage:
 * argv[1],[2],[3] matrix block size m, l, n
 * argv[4] OpenMP #Thread per task -- default 1
 * argv[5] the value of alpha -- default 1
 * argv[6] the value of beta -- default 0
 * */
int main(int argc, char* argv[]){
    if(argc!=7){
        std::cout<<"please check your input"<<std::endl;
    }
    int num_openmp_threads = atoi(argv[4]);
    double alpha = atof(argv[5]);
    double beta = atof(argv[6]);
    MPI_Init(nullptr,nullptr);
    int my_rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    int n_ranks = 0;
    MPI_Comm_size(MPI_COMM_WORLD,&n_ranks);
    // check if n_rank = k*k
    if((int)sqrt(n_ranks) * (int)sqrt(n_ranks) != n_ranks){
        std::cout<<"The number of MPI Process should be the square of a number."<<std::endl;
        MPI_Finalize();
        exit(-1);
    }
    int k = (int) sqrt(n_ranks);
    int m_size= atoi(argv[1]);
    int l_size= atoi(argv[2]);
    int n_size= atoi(argv[3]);
    if(my_rank==0){
        //generate three matrix, row-major
        std::vector<double> A(m_size*l_size*n_ranks,1.0);
        std::vector<double> B(l_size*n_size*n_ranks,1.0);
        std::vector<double> C(m_size*n_size*n_ranks,1.0);
    }

    //create the buffer
    std::vector<int> sub_A(m_size*l_size);
    std::vector<int> sub_B(l_size*n_size);
    std::vector<int> sub_C(m_size*n_size);
    if(my_rank==0){
        //re-order the data for the following MPI_Scatter
        std::vector<double> A_v(m_size*l_size);
        std::vector<double> B_v(l_size*n_size);
        std::vector<double> C_v(m_size*n_size);
        int A_count = 0, B_count = 0, C_count = 0;
        for (int p = 0; p < n_ranks; ++p) {
            //assign every sub-matrix, in my program, because all nodes except root node do not have sub-matrix
            //I just sent them the right sub-matrix after first shift. So they do not need to shift 1-k steps.
            int A_start = l_size *
            int B_start = n_size * l_start[p/k] + l_spilt[p/k]*n_start[p%k];
            for(int i = 0; i < m_spilt[p%k];++i){
                //copy one row
                std::copy(A.begin()+(m_start[p%]))
            }

        }

    }
    // other processes wait until root read all three matrix
    MPI_Barrier(MPI_COMM_WORLD);//I do not know if this is needed.


    //2. Permute Sub-Matrix
    MPI_Scatterv(B_v);
    int nbrs[4],dims[2]={k,k},periods[2]={1,1},reorder=0,coords[2];
    MPI_Comm cartcomm;
    MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,reorder,&cartcomm);
    MPI_Cart_coords(cartcomm,my_rank,2,coords);
    MPI_Cart_shift(cartcomm,0,1,&nbrs[0],&nbrs[1]);
    MPI_Cart_shift(cartcomm,1,1,&nbrs[2],&nbrs[3]);
    MPI_Sendrecv_replace(sub_A.data(),3232,MPI_DOUBLE,nbrs[2],3,nbrs[3],3,MPI_COMM_WORLD,&status);
    MPI_Sendrecv_replace(sub_B.data(),1323,MPI_DOUBLE,nbrs[0],4,nbrs[1],4,MPI_COMM_WORLD,&status);

    //3. Do k-times multiplication and movement, use OpenMP at the multiplication part.
    // multiplication
    for (int p = 0; p < k; ++p) {
        for (int j = 0; j < n_size; ++j) {
            for (int l = 0; l < l_size; ++l) {
                for (int i = 0; i < m_size; ++i) {
                    sub_C[i+j*m_size] += alpha * sub_A[i*l_size+l] * sub_B[l*n_size+j];
                }
            }
        }
    }
    //

    // I can not use MPI_sendrecv_replace if the sub-matrix size is different, however, I can achieve this by padding
}