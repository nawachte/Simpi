#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include "solveSystem.cpp"
#include <string.h>

void solve_REF_equations(double *M,double *answer,int dim){
  // int dim = M->get_x();
  int lastrc = dim-1;
  if (M[lastrc+dim*lastrc]!=0){
    answer[lastrc] = 0;
  }
  else{
    answer[lastrc] = 1;
  }
  for (int row=lastrc-1;row>=0;row--){
    int new_col = row;
    double known_sum = 0;
    for (int i=new_col+1;i<dim;i++){
      known_sum += (M[i+row*dim]*answer[i]);
    }
    if (known_sum==0 || M[new_col+row*dim]==0){
      answer[new_col] = 0;
    }
    else{
      // if (M[new_col+row*dim]==6){
      //   printf("6 known_sum: %.03f\n",known_sum);
      // }
      answer[new_col] = (-1*known_sum)/M[new_col+row*dim];
    }
  }
}

void swapRows(double *M,int i, int j,int dim){
  // int dim = M->get_x();
  for (int k=0; k<dim; k++){
      double temp = M[k+i*dim];
      M[k+i*dim] = M[k+j*dim];
      M[k+j*dim] = temp;
  }
}

void rowEchlon(double *M,int dim){
  // int dim = M->get_x();
  for(int k=0;k<dim;k++){
    int i_max = k;
    int v_max = M[k+i_max*dim];
    for (int i=k+1;i<dim;i++){
      if (abs(M[k+i*dim])>v_max){
        v_max = M[k+i*dim];
        i_max = i;
      }
    }
    if (i_max!=k){
      swapRows(M,k,i_max,dim);
    }
    for (int i=k+1;i<dim;i++){
      double f = M[k+i*dim]/M[k+k*dim];
      for (int j=k+1;j<dim;j++){
        M[j+i*dim] -= M[j+k*dim]*f;
        M[k+i*dim] = 0;
      }
    }
  }
  for (int i=0;i<dim*dim;i++){
    if (abs(M[i])<0.00001){
      M[i] = 0;
    }
  }
}

void eigenvector(matrix A, double *V, int par_id, int par_count){
  // printf("pi %d inside\n",par_id);
  // double *V;
   if (A.get_x()!=A.get_y()){
      throw std::invalid_argument("matrix must be a square matrix");
   }
  int dim = A.get_x();
  matrix *copyA = new matrix(dim,dim);
  for (int i=0;i<dim*dim;i++){
    copyA->arr[i] = A.arr[i];
  }
  double *eigval = eigenvalue(*copyA,par_id,par_count);

  // double *V;
  // int vectFD;
  // printf("%dgo to synch\n",par_id);
  // main_simpi->synch();
  // SIMPI_SYNCH();
  // printf("processes synched\n");
  // printf("id%d evals: %.0f %.0f %.0f\n",par_id,eigval[0],eigval[1],eigval[2]);



   // if (dim>par_id){
   if (par_id==0){
   //   printf("init shared mem\n");
   // vectFD = shm_open("eig_vect",O_RDWR|O_CREAT,0777);
   // if (vectFD==-1){
   //   printf("shared memory failed in eigen vector calculation\n");
   // }
   // ftruncate(vectFD,sizeof(double)*(dim*dim));
   // V = (double*)mmap(NULL,sizeof(double)*(dim*dim),PROT_READ|PROT_WRITE,MAP_SHARED,vectFD,0);
   // printf("shared mem started\n");
   // printf("do work process %d\n",par_id);
   // int work = dim/par_count;
   // if (work<1){
   //   work = 1;
   // }
   // int start = par_id*work;
   // int end = start+work;
   // int xtra = dim-(work*par_count);
   // if (xtra!=0 && par_count-xtra-1<par_id){
   //    start+=(par_id-(par_count-xtra));
   //    end = start+work+1;
   // }
   // for (int i=start;i<end;i++){
   // printf("start calculations\n");
   for (int i=0;i<dim;i++){
     // printf("\nvector calculations%d\n",i);
     // printf("init 0 vector");
     // vector *zero_vector = new vector(dim);
     // printf("declare double arr AminL\n");
     double AminL[dim*dim];
     // printf("double arr AminL declared\n");

     // double AminL[N][N+1];
     // printf("declare AminL\n");
     // matrix AminL = matrix(dim,dim);
     // printf("lambda = %.02f\n",eigval[i]);
     // printf("init AminL\n");
     for (int j=0;j<dim;j++){
       // zero_vector->arr[j] = 0;
       for (int k=0;k<dim;k++){
         if (j==k){
           AminL[k+dim*j] = A.arr[k+dim*j]-eigval[i];
         }
         else{
           AminL[k+j*dim] = A.arr[k+j*dim];
         }
       }
     }
     // printf("put in REF\n");
     rowEchlon(AminL,dim);
     // printf("out of REF\n");
     // vector *evect = new vector(dim);
     double evect[dim];
     // printf("solve REF equations\n");
     solve_REF_equations(AminL,evect,dim);
     for (int j=0; j<dim;j++){
       // printf("attempting to access shared mem\n");
       V[j+i*dim] = evect[j];
     }
   }
   }
   // printf("%d calculations finished\n",par_id);
   // // if (par_id!=0){
   // //   sleep(5);
   // //   printf("%d done sleep\n",par_id);
   // // }
   // if (par_id==0){
   // printf("inside vector method:\n");
   // std::cout<<"V address3: "<<V<<std::endl;
   // std::cout<<"|"<<V[0]<<" "<<V[1]<<" "<<V[2]<<"|"<<std::endl;
   // std::cout<<"|"<<V[3]<<" "<<V[4]<<" "<<V[5]<<"|"<<std::endl;
   // std::cout<<"|"<<V[6]<<" "<<V[7]<<" "<<V[8]<<"|"<<std::endl<<std::endl;
   // }
   // printf("id %d returning vectors\n",par_id);
}
