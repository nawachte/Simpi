#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

void pdDim_1_2(double *A, int dim, double *poly, int *lambdas){
   /*
 dim == 1 only when polyDet is called from the eigen value method,
 this will most likely not exist once the eigen value method is made and
 1x1 matrices are easily calculated on their own.
   */
 if (dim == 1){
     poly[0] = A[0];
     poly[1] = -1;
     return;
 }
 // base case
else {
     double polyad[3];
     double polybc[3];
     //mulitply ad
     if (lambdas[0] == 0 && lambdas[3]==0){
          polyad[0] = A[0]*A[3];
          polyad[1] = 0;
          polyad[2] = 0;
     }
     else if(lambdas[0] == 1 && lambdas[3]==0){
          polyad[0] = A[0]*A[3];
          if (A[3]==0){
             polyad[1] = 0;
          }
          else{
             polyad[1] = -1*A[3];
          }
          polyad[2] = 0;
     }
     else if(lambdas[0] == 0 && lambdas[3]==1){
          polyad[0] = A[0]*A[3];
          if (A[0] == 0){
             polyad[1] = 0;
          }
          else{
             polyad[1] = -1*A[0];
          }
          polyad[2] = 0;
     }
     else{
          polyad[0] = A[0]*A[3];
          if (A[0]==0 && A[3]==0){
             polyad[1] = 0;
          }
          else{
             polyad[1] = -1*(A[0]+A[3]);
          }
          polyad[2] = 1;
     }
     //multiply bc
     if (lambdas[1] == 0 && lambdas[2]==0){
          polybc[0] = A[1]*A[2];
          polybc[1] = 0;
          polybc[2] = 0;
     }
     else if(lambdas[1] == 1 && lambdas[2]==0){
          polybc[0] = A[1]*A[2];
          if (A[2]==0){
            polybc[1] = 0;
          }
          else{
            polybc[1] = -1*A[2];
          }
          polybc[2] = 0;
     }
     else if(lambdas[1] == 0 && lambdas[2]==1){
          polybc[0] = A[1]*A[2];
          if (A[1]==0){
             polybc[1] = 0;
          }
          else{
             polybc[1] = -1*A[1];
          }
          polybc[2] = 0;
     }
     else{
          polybc[0] = A[1]*A[2];
          if (A[1]==0 && A[2]==0){
             poly[1] = 0;
          }
          else{
             polybc[1] = -1*(A[1]+A[2]);
          }
          polybc[2] = 1;
     }
     // subtracting ad-bc and putting result in poly
     if (polybc[0] == 0){
        poly[0] = polyad[0];
     }
     else{
        poly[0] = polyad[0]-polybc[0];
     }
     if (polybc[1]==0){
        poly[1] = polyad[1];
     }
     else{
        poly[1] = polyad[1]-polybc[1];
     }
     if (polybc[2]==0){
        poly[2] = polyad[2];
     }
     else{
        poly[2] = polyad[2]-polybc[2];
     }
     return;
 }
}

void polyDet_aux(double *A, int dim, double *poly, int *lambdas){
  if (dim == 1 || dim == 2){
     pdDim_1_2(A,dim,poly,lambdas);
     return;
  }
   int sign = 0; //0 to add 1 to subtract
   // for each value in the top row
   for (int x=0;x<dim;x++){
      // make new matrix
      double subA[(dim-1)*(dim-1)];
      int idx = 0;
      int newLambdas[(dim-1)*(dim-1)];
      for (int i=1;i<dim;i++){
           for (int j=0;j<dim;j++){
               if (j!=x){
                   subA[idx] = A[j+dim*i];
                   if (lambdas[j+dim*i]==1){
                       newLambdas[idx] = 1;
                   }
                   else{
                       newLambdas[idx] = 0;
                   }
                   idx++;
               }
           }
      }
      double subPoly[dim+1];
      for (int sp=0;sp<dim+1;sp++){
          subPoly[sp] = 0;
      }
      polyDet_aux(subA,dim-1,subPoly,newLambdas);
      // add and subtract sub-determinants
      if (sign == 0){
           for (int i=0;i<sizeof(subPoly)/sizeof(double);i++){
               poly[i] += subPoly[i]*A[x];
           }
           if (lambdas[x]==1){
               for (int i=0;i<sizeof(subPoly)/sizeof(double);i++){
                   poly[i+1] -= subPoly[i];
               }
           }
           else{
               // highest possible power is of lambda is len(poly)-2
               poly[dim+1] = 0;
           }
           sign = 1;
      }
      else{
           double subtractPoly[dim+1];
           for (int i=0;i<dim+1;i++){
               subtractPoly[i] = 0;
           }
           for (int i=0;i<sizeof(subPoly)/sizeof(double);i++){
               subtractPoly[i] += subPoly[i]*A[x];
           }
           if (lambdas[x]==1){
               for (int i=0;i<sizeof(subPoly)/sizeof(double);i++){
                   subtractPoly[i+1] -= subPoly[i];
               }
           }
           else{
               // highest possible power is of lambda is len(poly)-2
               subtractPoly[dim+1] = 0;
           }
           for (int i=0;i<(sizeof(subPoly)/sizeof(double))+1;i++){
               poly[i] -= subtractPoly[i];
           }
           sign = 0;
      }
   }
}

void pre_polyDet_aux(double *A, int dim, double *poly, int *lambdas, int sign, int x){
  double subA[(dim-1)*(dim-1)];
  int idx = 0;
  int newLambdas[(dim-1)*(dim-1)];
  for (int i=1;i<dim;i++){
      for (int j=0;j<dim;j++){
          if (j!=x){
              subA[idx] = A[j+dim*i];
              if (lambdas[i+dim*j]==1){
                  newLambdas[idx] = 1;
              }
              else{
                  newLambdas[idx] = 0;
              }
              idx++;
          }
      }
  }
  double subPoly[dim+1];
  for (int sp=0;sp<dim+1;sp++){
     subPoly[sp] = 0;
  }
  polyDet_aux(subA,dim-1,subPoly,newLambdas);
  // add and subtract sub-determinants
  if (sign == 0){
      for (int i=0;i<sizeof(subPoly)/sizeof(double);i++){
          poly[i] += subPoly[i]*A[x];
      }
      if (lambdas[x]==1){
          for (int i=0;i<sizeof(subPoly)/sizeof(double);i++){
              poly[i+1] -= subPoly[i];
          }
      }
      else{
          // highest possible power is of lambda is len(poly)-2
          poly[dim+1] = 0;
      }
  }
  else{
      double subtractPoly[dim+1];
      for (int i=0;i<dim+1;i++){
          subtractPoly[i] = 0;
      }
      for (int i=0;i<sizeof(subPoly)/sizeof(double);i++){
          subtractPoly[i] += subPoly[i]*A[x];
      }
      if (lambdas[x]==1){
          for (int i=0;i<sizeof(subPoly)/sizeof(double);i++){
              subtractPoly[i+1] -= subPoly[i];
          }
      }
      else{
          // highest possible power is of lambda is len(poly)-2
          subtractPoly[dim+1] = 0;
      }
      for (int i=0;i<(sizeof(subPoly)/sizeof(double))+1;i++){
          poly[i] -= subtractPoly[i];
      }
  }
}

double power(double num, int pow){
   double temp = num;
   for (int i=0;i<pow-1;i++){
      temp=temp*num;
   }
   return temp;
}
/*
The input poly is a float array with each index representing a multiple of that
index-power of x in the polynomial equation d+cx+bx^2+..+ax^(len(P)-1)

The method polyRoots solve for the roots of this equation and returns a float
array of all zeroes of P

This method is based off of The Rational Zero Theorem
*/
void polyRoots(double *poly,double *roots,int polySize){
    /*
      number roots cant exceed number of factors of the constant term of the polynomial which is always less than
      or equal to the number itself
      the rest will be filled in with 0s as 0 cannot be a valid eigan value

      P = factors of the constant term
      Q = factors of the leading coefficient
    */
    // round decimal place
    for (int i=0;i<polySize;i++){
      poly[i] = floor(poly[i]*10000)/10000;
    }
    if (poly[0]==0){
      // take out a factor of lambda and find the new solution with a root of zero
      double subpoly[polySize-1];
      for (int i=1;i<polySize;i++){
         subpoly[i-1] = poly[i];
      }
      double subroots[polySize-1];
      polyRoots(subpoly,subroots,polySize-1);
      roots[0] = 0;
      for (int i=1;i<polySize;i++){
         roots[i] = subroots[i-1];
      }
      return;
    }
    int P = (int)abs(poly[0]);
    int Q = (int)abs(poly[polySize-1]);
    // decimal adjustment
    int factor = 0;
    while (abs(P-abs(poly[0]))>0.00001 || abs(Q-abs(poly[polySize-1]))>0.00001){
      for (int i=0;i<polySize;i++){
         poly[i] = 10*poly[i];
      }
      factor++;
      P = (int)abs(poly[0]);
      Q = (int)abs(poly[polySize-1]);
   }
    /*
     these max's are causing a limit on how many roots can be found
     but they are causing an error when creating the factorsofX arrays
     when they are too big
    */
    int maxnumPs = P/2;
    int maxnumQs = Q/2;
    // find factors of P and Q
    int factorsOfP[maxnumPs];
    int fopIDX = 0;
    for (int i=1;i<P+1;i++){
        if (P%i == 0){
            factorsOfP[fopIDX] = i;
            fopIDX++;
        }
    }
    int factorsOfQ[maxnumQs];
    int foqIDX = 0;
    for (int i=1;i<Q+1;i++){
        if (Q%i == 0){
            factorsOfQ[foqIDX] = i;
            foqIDX++;
        }
    }
    // create our P/Q values
    int pos_rootsIDX = 0;
    double possible_roots[foqIDX*fopIDX*2];
    for (int i=0;i<foqIDX;i++){
        for (int j=0;j<fopIDX;j++){
            possible_roots[pos_rootsIDX] = (double)factorsOfP[j]/factorsOfQ[i];
            pos_rootsIDX++;
            possible_roots[pos_rootsIDX] = -1*(double)factorsOfP[j]/factorsOfQ[i];
            pos_rootsIDX++;
        }
    }
    // restore poly if necessary---testing
    if (factor!=0){
      for (int i=0;i<factor;i++){
         for (int j=0;j<polySize;j++){
            poly[j] = poly[j]/10;
         }
      }
   }
    // evaluate possible roots
    int rootsIDX = 0;
    for (int i=0;i<pos_rootsIDX;i++){
        double evaluation;
        for (int j=0;j<polySize;j++){
            if (j==0){
                evaluation = poly[0];
            }
            else{
                evaluation += (double)poly[j]*power(possible_roots[i],j);
            }
        }
        int already_root = 0;
        for (int j=0;j<rootsIDX;j++){
           if (roots[j]==possible_roots[i]){
             already_root = 1;
          }
        }
        if (abs(evaluation)<0.00001 && already_root==0){

            roots[rootsIDX] = possible_roots[i];
            rootsIDX++;
        }
    }
}

void polyDet(int polyFD, int par_id, int par_count, int dim, double *A, double *poly, int *lambdas){
   if (dim<3){
      pdDim_1_2(A,dim,poly,lambdas);
      return;
   }
   int work = dim/par_count;
   if (work<1){work=1;}
   int start = par_id*work;
   int end = start+work;
   int xtra = dim-(work*par_count);
   if (xtra!=0 && par_count-xtra-1<par_id){
      start+=(par_id-(par_count-xtra));
      end = start+work+1;
   }
   // split up the remaining processes
   for (int x=start;x<end;x++){
      int sign = 0;
      if (x%2!=0){
         sign = 1;
      }
      pre_polyDet_aux(A,dim,poly,lambdas,sign,x);
   }
}
double *eigenvalue(matrix Anaught,int par_id, int par_count){
   // return Anaught.arr;
   int dim = Anaught.get_x();
   if (Anaught.get_x()!=Anaught.get_y()){
      throw std::invalid_argument("matrix must be a square matrix\n");
   }
   // create lambda matrix and copy of A
   int lambdas[dim*dim];
   matrix A = matrix(dim,dim);
   for (int i=0;i<dim*dim;i++){
      A.arr[i] = Anaught.arr[i];
      if (i%(dim+1)){
         lambdas[i] = 0;
      }
      else{
         lambdas[i] = 1;
      }
   }
   // init poly memory
   double *poly;
   double *roots;
   int polyFD;
   int rootsFD;
   if (par_id==0){
      polyFD = shm_open("eig_val_poly", O_RDWR|O_CREAT, 0777);
      if (polyFD==-1){
         printf("shared memory failed in eigen value calculation\n");
      }
      ftruncate(polyFD,sizeof(double)*(dim+1));
      poly = (double*)mmap(NULL,sizeof(double)*(dim+1),PROT_READ|PROT_WRITE,MAP_SHARED,polyFD,0);
      rootsFD = shm_open("eig_val_roots", O_RDWR|O_CREAT, 0777);
      if (rootsFD==-1){
         printf("shared memory failed in eigen value calculation\n");
      }
      ftruncate(rootsFD,sizeof(double)*(dim+1));
      roots = (double*)mmap(NULL,sizeof(double)*(dim+1),PROT_READ|PROT_WRITE,MAP_SHARED,rootsFD,0);
      for (int i=0;i<dim+1;i++){
         poly[i] = 0;
         roots[i] = 0;
      }
   }
   main_simpi->synch();
   if (par_id!=0){
      // sleep(2);
      polyFD = shm_open("eig_val_poly", O_RDWR, 0777);
      poly = (double*)mmap(NULL,sizeof(double)*(dim+1),PROT_READ|PROT_WRITE,MAP_SHARED,polyFD,0);
      rootsFD = shm_open("eig_val_roots", O_RDWR, 0777);
      roots = (double*)mmap(NULL,sizeof(double)*(dim+1),PROT_READ|PROT_WRITE,MAP_SHARED,rootsFD,0);
   }
   // get the polynomial from the determinent
   if (dim>par_id){
   polyDet(polyFD, par_id, par_count, dim, A.arr, poly, lambdas);
   }
   main_simpi->synch();
   // double *roots = (double*)mmap(NULL,sizeof(double)*dim+1,PROT_READ|PROT_WRITE,MAP_SHARED|MAP_ANONYMOUS,-1,0);
   // for (int i=0; i<dim+1; i++){
   //    roots[i] = 0;
   // }
   if (par_id==0){
      polyRoots(poly,roots,dim+1);
   }
   shm_unlink("eig_val_poly");
   // if (par_id==0){
   // printf("A2\n");
   //   for (int i=0; i<dim*dim;i++){
   //       printf("%.2f ",Anaught.arr[i]);
   //       if ((i+1)%dim==0 ){
   //        printf("\n");
   //       }
   //   }
   // }
   main_simpi->synch();
   return roots;
}
