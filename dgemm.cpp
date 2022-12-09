//
// Created by Hao Feng on 11/28/22.
//

#include "iostream"
#include <omp.h>
#include "mpi.h"
#include <chrono>                   // std::chrono::high_resolution_clock::now()
#include <cmath>                    // std::sin, M_PI
#include <numeric>                  // std::transform_reduce
#include <vector>

/*
 * args usage:
 * argv[1],[2],[3] matrix block size m, l, n
 * argv[4] the value of alpha -- default 1
 * argv[5] the value of beta -- default 0
 * */
int main(int argc, char* argv[]){
    if(argc!=6){
        std::cout<<"please check your input"<<std::endl;
    }
    double alpha = atof(argv[4]);
    double beta = atof(argv[5]);
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
    int mm = m_size*k;
    int nn = n_size*k;
    int ll = l_size*k;
    //create the buffer
    std::vector<int> sub_A(m_size*l_size);
    std::vector<int> sub_B(l_size*n_size);
    std::vector<int> sub_C(m_size*n_size);
    std::vector<double> A_v(0);
    std::vector<double> B_v(0);
    std::vector<double> C_v(0);
    if(my_rank==0){
        //generate three matrix, row-major
        std::cout<<"root setup start"<<std::endl;
        std::vector<double> A(m_size*l_size*n_ranks,1.0);
        std::vector<double> B(l_size*n_size*n_ranks,1.0);
        std::vector<double> C(m_size*n_size*n_ranks,1.0);
        //re-order the data for the following MPI_Scatter
        A_v.resize(m_size*l_size*n_ranks);
        B_v.resize(l_size*n_size*n_ranks);
        C_v.resize(m_size*n_size*n_ranks);
        int A_count = 0, B_count = 0, C_count = 0;
        for (int p = 0; p < n_ranks; ++p) {
            //assign every sub-matrix, in my program, because all nodes except root node do not have sub-matrix
            //I just sent them the right sub-matrix after first shift. So they do not need to shift 1-k steps.
            int ii = p/k;
            int jj = p%k;
            int A_start = ii*k*(m_size*l_size) + jj*l_size;
            int B_start = ii*k*(l_size*n_size) + jj*n_size;
            int C_start = ii*k*(m_size*n_size) + jj*n_size;
            for(int i = 0; i < m_size;++i){
                //copy one row
                std::copy(A.begin()+A_start+i*ll,A.begin()+A_start+i*ll+l_size,A_v.begin()+A_count);
                A_count += l_size;
            }
            for (int i = 0; i < l_size; ++i) {
                std::copy(B.begin()+B_start+i*nn,B.begin()+B_start+i*nn+n_size,B_v.begin()+B_count);
                B_count += n_size;
            }
            for (int i = 0; i < m_size; ++i) {
                std::copy(C.begin()+C_start+i*nn,C.begin()+C_start+i*nn+n_size,C_v.begin()+C_count);
                C_count += n_size;
            }
        }
        std::cout<<"root setup finished"<<std::endl;
    }

    //2. Permute Sub-Matrix
    // other processes wait until root read all three matrix
    MPI_Barrier(MPI_COMM_WORLD);//I do not know if this is needed.
    std::cout<<"can reach there\n";
    MPI_Scatter(A_v.data(),m_size*l_size,MPI_DOUBLE,sub_A.data(),m_size*l_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatter(B_v.data(),n_size*l_size,MPI_DOUBLE,sub_B.data(),n_size*l_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatter(C_v.data(),m_size*n_size,MPI_DOUBLE,sub_C.data(),m_size*n_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
    std::cout<<"rank"<<my_rank<<"can reach there\n";
    /*
    int nbrs[4],dims[2]={k,k},periods[2]={1,1},reorder=0,coords[2];
    MPI_Comm cartcomm;
    MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,reorder,&cartcomm);
    MPI_Cart_coords(cartcomm,my_rank,2,coords);
    MPI_Cart_shift(cartcomm,0,1,&nbrs[0],&nbrs[1]);
    MPI_Cart_shift(cartcomm,1,1,&nbrs[2],&nbrs[3]);
    //MPI_Sendrecv_replace(sub_A.data(),m_size*l_size,MPI_DOUBLE,nbrs[2],3,nbrs[3],3,MPI_COMM_WORLD,&status);
    //MPI_Sendrecv_replace(sub_B.data(),l_size*n_size,MPI_DOUBLE,nbrs[0],4,nbrs[1],4,MPI_COMM_WORLD,&status);
    */
    //3. Do k-times multiplication and movement, use OpenMP at the multiplication part.
    // multiplication

#pragma omp parallel for
    for (int i = 0; i < m_size*n_size; ++i) {
        sub_C[i] *= beta;
    }
    /*
    for (int p = 0; p < k; ++p) {
#pragma omp parallel for
        for (int j = 0; j < n_size; ++j) {
            for (int l = 0; l < l_size; ++l) {
                for (int i = 0; i < m_size; ++i) {
                    sub_C[i+j*m_size] += alpha * sub_A[i*l_size+l] * sub_B[l*n_size+j];
                }
            }
        }
    }
     */
    // Gather all the result
    MPI_Finalize();
}