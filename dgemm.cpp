//
// Created by Hao Feng on 11/28/22.
//
/*
 * args usage:
 * argv[1] matrix file name
 * argv[2] MPI world size -- must can be squared
 * argv[3] OpenMP #Thread per task -- default 1
 * */
int main(int argc, char* argv[]){
    int k = atoi(argv[2]);
    //1. Reorder Matrix and setup MPI


    //2. Permute Sub-Matrix

    //3. Do k-times multiplication and movement, use OpenMP at the multiplication part.
}