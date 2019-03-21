#include <mex.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include "shape.h"
//#include <omp.h>


void fluxHorizontal(int M, int N, double *u, double *v, double *psi ,double *nx,
        double *ny, double *alpha, double dx, double dt, double *psi2){
    int i,j,k,m,n, shx, shy, co = 0; double Lx, Ly, Lm, L;
    for(i = 0; i < (M+4)*(N+4); i++)
        psi2[i] = 0.0;
    double X0, X1, Y0, Y1;
    for(i = 2; i < M+2; i++){//for(i = M + 1; i > 1; i--)
        for(j = 2; j < N + 2; j++){
            int condition1 = (alpha[i+(M+4)*(j)]>1.5e19)&&(alpha[i+(M+4)*(j)]<2.5e19);
            int condition2 = psi[i+(M+4)*(j)]>1e-11;
            int condition3 = 0;
            for (shx = -2; shx < 3; shx++){ 
                for (shy = -2; shy < 3; shy++){
                    condition3 = condition3||(alpha[i+shx+(M+4)*(j+shy)]<1e19);
                }
            }
            if (condition3){
                double ue = -u[i-1 + (M+1)*(j-1)]*dt/dx, uw = -u[i-2 + (M+1)*(j-1)]*dt/dx, un = -v[i-1 + (M+2)*(j-1)]*dt/dx, us = -v[i-1 + (M+2)*(j-2)]*dt/dx,
                        x0 = (u[i-2 + (M+1)*(j-1)]+u[i-2 + (M+1)*(j-2)])*dt/(2.0*dx)-.5, y0 = (v[i-1 + (M+2)*(j-2)]+v[i-2 + (M+2)*(j-2)])*dt/(2.0*dx)-.5,
                        x2 = (u[i-1 + (M+1)*(j-1)]+u[i-1 + (M+1)*(j-2)])*dt/(2.0*dx)+.5, y2 = (v[i   + (M+2)*(j-2)]+v[i-1 + (M+2)*(j-2)])*dt/(2.0*dx)-.5,
                        x4 = (u[i-1 + (M+1)*(j  )]+u[i-1 + (M+1)*(j-1)])*dt/(2.0*dx)+.5, y4 = (v[i   + (M+2)*(j-1)]+v[i-1 + (M+2)*(j-1)])*dt/(2.0*dx)+.5,
                        x6 = (u[i-2 + (M+1)*(j  )]+u[i-2 + (M+1)*(j-1)])*dt/(2.0*dx)-.5, y6 = (v[i-1 + (M+2)*(j-1)]+v[i-2 + (M+2)*(j-1)])*dt/(2.0*dx)+.5;
                double a[] = {0, x0 - x2, 0, x2 - x4, 0, x4 - x6, 0, x6 - x0};
                double b[] = {y2/2 - y0/2,y0-y2,y2*(1.0/2.0)-y4*(1.0/2.0),y2-y4,y4*(1.0/2.0)-y6*(1.0/2.0),y4-y6,y0*(1.0/2.0)-y6*(1.0/2.0), -y0+y6};
                double c[] = {x0/2 - x2/2,0, x4/2 - x2/2, 0, x6/2 - x4/2, 0, x6/2 - x0/2};
                double d[] = {us - x0/4 + x2/4 + y0/4 + y2/4 + 1.0/4,
                x0*x0/2 - x2*x2/2 + y0*y0/2 - y2*y2/2,
                ue + x2/4 + x4/4 + y2/4 - y4/4 - 1.0/4,
                x2*x2/2 - x4*x4/2 + y2*y2/2 - y4*y4/2,
                un - x4/4 + x6/4 + y4/4 + y6/4 - 1.0/4,
                x4*x4/2 - x6*x6/2 + y4*y4/2 - y6*y6/2,
                uw + x0/4 + x6/4 - y0/4 + y6/4 + 1.0/4,
                - x0*x0/2 + x6*x6/2 - y0*y0/2 + y6*y6/2};
                
                for(k = 0;k < 8;k+=2) solve(a,b,c,d,k,k+2);
                
                double x1 = d[0], y1 = d[1], x3 = d[2], y3 = d[3], x5 = d[4], y5 = d[5], x7 = d[6], y7 = d[7];
                if (j==2) {x1 = 0; y1 = -.5;}
                //if (fabs(ue) < 1e-12) {x3 = .5; y3 = 0;}
                if (j==N+1) {x5 = 0; y5 = .5;}
                //if (fabs(uw) < 1e-12) {x7 = -.5; y7 = 0;}
                
                double X[] = {x0, x1, x2, x3, x4, x5, x6, x7};
                double Y[] = {y0, y1, y2, y3, y4, y5, y6, y7};
                double area = 1/2*(x0*y1+x1*y2+x2*y3+x3*y4+x4*y5+x5*y6+x6*y7+x7*y0-
                        (y0*x1+y1*x2+y2*x3+y3*x4+y4*x5+y5*x6+y6*x7+y7*x0));
                //psi2[i+(M+4)*(j)] = 0;
                
                /*for(k = 0; k < 8; k++){
                    preview[co+(M*N)*(k)] = X[k]+i;
                    preview[co+(M*N)*(k+8)] = Y[k]+j;
                }
                co++;*/
                
                double psit[3][3], mx[3][3], my[3][3], alpha1[3][3]; int in;
                poly p1[3][3], pmain = poly_new(), res;  vec_t tmp;
                for (k = 0; k < 8; k++){
                    vnew(X[k], Y[k], &tmp);
                    poly_append(pmain, &tmp);
                }
                
                for(m = 0; m < 3;m++){
                    for( n = 0; n < 3; n++){
                        in = i+m-1+(M+4)*(j+n-1);
                        p1[m][n] = poly_new();
                        psit[m][n] = psi[in];
                        mx[m][n] = nx[in];
                        my[m][n] = ny[in];
                        alpha1[m][n] = alpha[in];
                    }
                }
                
                for(m = -1; m < 2; m++){
                    for(n = -1; n < 2; n++){
                        if ((alpha1[m+1][n+1] < 1e19)||(psit[m+1][n+1] > 0.5)){
                            subComputePolygon(mx[m+1][n+1],my[m+1][n+1],alpha1[m+1][n+1],p1[m+1][n+1],m,n);
                            //
                            //poly_t clipper = {sizeof(p1[m+1][n+1]->v)/sizeof(vec_t), 0, p1[m+1][n+1]->v};
                            //
                            res = poly_clip(pmain, p1[m+1][n+1]);
                            psi2[i+(M+4)*(j)] += polygonArea(res);
                            poly_free(res);
                        }
                        poly_free(p1[m+1][n+1]);
                    }
                }
                poly_free(pmain);
                ///end/////
            }
            //else if (condition1) psi2[i+(M+4)*(j)] = 1.0;
            //else psi2[i+(M+4)*(j)] = 0.0;
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    srand(554);
    /* Macros for the ouput and input arguments */
    double *B, *U, *V, *mxin,  *myin, *a1in, *psin, dx, dy, dt;
    int M, N, i, j, k, num;
    M = mxGetScalar(prhs[0]); /* Get the dimensions of A */
    N = mxGetScalar(prhs[1]);
    U    = mxGetPr(prhs[2]);
    V    = mxGetPr(prhs[3]);
    mxin = mxGetPr(prhs[4]);
    myin = mxGetPr(prhs[5]);
    a1in = mxGetPr(prhs[6]);
    psin = mxGetPr(prhs[7]);
    dx = mxGetScalar(prhs[8]); dy = dx;
    dt = mxGetScalar(prhs[9]);
    
    plhs[0] = mxCreateDoubleMatrix(M+4, N+4, mxREAL); /* Create the output matrix */
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    
    B  = mxGetPr(plhs[0]);
    double *u = malloc(sizeof(double)*(M+1)*(N+2));
    double *v = malloc(sizeof(double)*(M+2)*(N+1));
    for (i = 0; i < (M+1)*(N+2); i++) u[i] = -U[i];
    for (i = 0; i < (M+2)*(N+1); i++) v[i] = -V[i];
    /***********************************************************
     *Begin Computation
     ***********************************************************/
    
    
    clock_t start, end;
    double* cpu_time_used = mxGetPr(plhs[1]);
    start = clock();
    
    fluxHorizontal(M,N,u,v,psin,mxin,myin,a1in,dx,dt,B);
    
    end = clock();
    cpu_time_used[0] = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    free(u); free(v);
    /*nx and ny are not normal but sum vectors*/
    return;
}



