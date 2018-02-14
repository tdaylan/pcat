#include <stdio.h>
#include <stdlib.h>
//#include "mkl_example.h"
//#include "mkl_cblas.h"
//#include "export.h"
//#include "i_malloc.h"

//void myprint(void);

void myprint(){
    printf("hello world\n");
}

//void retr_deflsubh(int argc, void * argv[]){
//
//    long nstar,n,nc,k;
//    int      i,i2,j,j2,rad,istar,xx,yy;
//    float    *A, *B, *C, alpha, beta;
//    float    (*image)[100][100];         // image size hardwired
//    int      (*x)[10000], (*y)[10000];   // max number of stars 10000
//
//    /* allocate pointers from IDL */
//    nstar = *((IDL_LONG *)argv[0]);
//    nc    = *((IDL_LONG *)argv[1]);
//    k     = *((IDL_LONG *)argv[2]);
//    A     = (float *)argv[3];            // npar by nstar matrix
//    B     = (float *)argv[4];            // npix by npar coeff matrix
//    C     = (float *)argv[5];
//    x     = argv[6];
//    y     = argv[7];
//    image = argv[8];
//
//    n = nc*nc;
//    rad = nc/2;
//
//    alpha = 1.0; beta = 0.0;
//
////  matrix multiplication
//    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nstar, n, k, alpha, A, k, B, n, beta, C, n);
//
////    printf("       M=%i  N=%i  K=%i\n", m, n, k);
//
////  loop over stars, insert psfs into image    
//    for (istar = 0 ; istar < nstar ; istar++)
//    {
//	xx = (*x)[istar];
//	yy = (*y)[istar];
//	for (j = yy-rad, j2 = istar*nc*nc ; j <= yy+rad ; j++, j2+=nc)
//	    for (i2 = 0, i = xx-rad ; i2 < nc ; i2++, i++)
//		(*image)[j][i] += C[i2+j2];
//    }
//
//}
//
