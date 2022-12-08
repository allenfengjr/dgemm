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
 * argv[1],[2],[3] matrix file name for A(m*l), B(l*n), C(m*n)
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
    int m_size, n_size, l_size;
    if(my_rank==0){
        MMIOReader A(argv[1]);
        MMIOReader B(argv[2]);
        MMIOReader C(argv[3]);
        //check shape of three matrix
        if(A.n_rows!=C.n_rows||B.n_columns!=C.n_columns||A.n_columns!=B.n_rows){
            std::cout<<"the shapes of matrix can not match"<<std::endl;
            MPI_Finalize();
            exit(-1);
        }
        m_size = A.n_rows;
        l_size = A.n_columns;
        n_size = B.n_columns;
        MPI_Bcast(&m_size,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&n_size,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&l_size,1,MPI_INT,0,MPI_COMM_WORLD);
    }
    std::vector<int> m_spilt(k);
    std::vector<int> n_spilt(k);
    std::vector<int> l_spilt(k);
    for (int i = 0; i < k; ++i) {
        m_spilt[i] = (i<(int)m_size%k)? ceil(m_size/k): floor(m_size/k);
        n_spilt[i] = (i<(int)n_size%k)? ceil(n_size/k): floor(n_size/k);
        l_spilt[i] = (i<(int)l_size%k)? ceil(l_size/k): floor(l_size/k);
    }
    //create the largest space, if
    std::vector<int> sub_A(m_spilt[0]*l_spilt[0]);
    std::vector<int> sub_B(l_spilt[0]*n_spilt[0]);
    std::vector<int> sub_C(m_spilt[0]*n_spilt[0]);
    if(my_rank==0){
        //re-order the data for the following MPI_Scatter
        std::vector<double> A_v(m_size*l_size);
        std::vector<double> B_v(l_size*n_size);
        std::vector<double> C_v(m_size*n_size);
        for (int i = 0; i < n_ranks; ++i) {
            //assign every sub-matrix, in my program, because all nodes except root node do not have sub-matrix
            //I just sent them the right sub-matrix after first shift. So they do not need to shift 1-k steps.
            

        }

    }
    // other processes wait until root read all three matrix
    MPI_Barrier(MPI_COMM_WORLD);//I do not know if this is needed.


    //2. Permute Sub-Matrix
    MPI_Scatterv(B_v);


    //3. Do k-times multiplication and movement, use OpenMP at the multiplication part.
    // multiplication
    for(int i = 0 ; i < ; ++i){
        for (int j = 0; j < ; ++j) {

        }
    }
    //

    // I can not use MPI_sendrecv_replace if the sub-matrix size is different, however, I can achieve this by padding
}