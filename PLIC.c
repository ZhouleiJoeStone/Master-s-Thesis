#include "mex.h"
#include "math.h"
#include "stdint.h"
#include <time.h>
#include <stdbool.h>
double tol = 1e-5;
double sign(double x) {return ((x > 0) ? 1.0 : ((x < 0) ? -1.0 : 0.0));}
double step(double x) {return ((x > 0) ? 1.0 : ((x < 0) ?  0.0 : 0.5));}
double min(double a, double b) {return ((a > b) ? b : a);}
double max(double a, double b) {return ((a > b) ? a : b);}
int caseno; //Choice of Reconstruction
//1: FDM, 2: COM, 3: ELVIRA, 4: LVIRA
double pg[3][3]; //global variable
double eps(double x);
double findmin(double x0);
double alphacomp(double psi, double mx1, double my1){     
    double a = 0, m1, m2, psi1, chi, snrm, mx, my; //computing alpha for mx
    snrm = fabs(mx1)+fabs(my1);
    mx = mx1/snrm;
    my = my1/snrm;
    if ((psi<1.0e-12) || (psi>(1.0-1.0e-12))){
        a = 1.0e20;
    }
    else{
        m1 = min(fabs(mx),fabs(my));
        m2 = max(fabs(mx),fabs(my));
        psi1 = m1/(2.0*(1.0-m1));
        chi = min(psi,1.0-psi);
        if (chi<psi1)
            a = sqrt(2.0*m1*(1-m1)*chi);
        if ((chi>=psi1) && (chi<=.5))
            a = (1-m1)*chi+m1/2.0;
        if (psi>.5)
            a = 1.0-a;
    }
    return a;
}

double subalphacomp(double psi, double mx1, double my1){
    psi = 1.0 - psi;     
    double a = 0, m1, m2, psi1, chi, snrm, mx, my; //computing alpha for mx
    snrm = fabs(mx1)+fabs(my1);
    mx = mx1/snrm;
    my = my1/snrm;
    if ((psi<0) || (psi>(1.0))){
        a = 1.0e20;
    }
    else{
        m1 = min(fabs(mx),fabs(my));
        m2 = max(fabs(mx),fabs(my));
        psi1 = m1/(2.0*(1.0-m1));
        chi = min(psi,1.0-psi);
        if (chi<psi1)
            a = sqrt(2.0*m1*(1-m1)*chi);
        if ((chi>=psi1) && (chi<=.5))
            a = (1-m1)*chi+m1/2.0;
        if (psi>.5)
            a = 1.0-a;
        a = a + min(0,mx) + min(0,my);
        a = (a-0.5*(mx+my))/sqrt(mx*mx+my*my);
    }
    return a;
}

double areacomp(double a, double b, double c){
    //Compute Area of PLIC
    if (a==0||b==0){
        a+=1e-10;
        b+=1e-10;
    }
    double capcb = c*a*a + c*b*b;
    double hp = (a+b)/2.0;
    double hm = (a-b)/2.0;
    
    double cpp = capcb + hp;
    double cmp = capcb - hp;
    double cpm = capcb + hm;
    double cmm = capcb - hm;
    
    double area = 0.5+(cpp*cpp*sign(cpp) + cmp*cmp*sign(cmp) - cpm*cpm*sign(cpm) - cmm*cmm*sign(cmm))/(4.0*a*b);
    return area;
}

double findError(double t){
  int i,j; double sum=0, ps[3][3];
  double cost = cos(t), sint = sin(t);
  double a = subalphacomp(pg[1][1], cost, sint);
  for(j = -1; j <= 1; j++){
    for(i = -1; i <= 1; i++){
      if (!((i==0)&&(j==0)))
      ps[i+1][j+1] = areacomp(cost, sint, -a+j*sint+i*cost);
      else ps[1][1] = pg[1][1];
    }
  }
  for (i = 0; i < 3; i++){
    for (j = 0; j < 3; j++){
      sum+=pow(pg[i][j]-ps[i][j],2.0);
    }
  }
  return sum;
}

void normalvec(int M, int N,  double dx, double dy, double p[M+4][N+4], 
        double nx[M+4][N+4], double ny[M+4][N+4], double alpha[M+4][N+4]){
    int16_t i,j,m,n;

    
    for (i = 0; i < M+4; i++){
        for(j = 0; j < N+4; j++){
            int BC = (i==0) || (i==M+3) || (j==0) || (j==N+3);
            int BC2= (i==1) || (i==M+2) || (j==1) || (j==N+2);
            
            
            if ((/*((mx[i][j]==0) && (my[i][j]==0)) ||*/ BC) && ((p[i][j]<1e-8) || (p[i][j]>(1-1e-8)))){
                nx[i][j] = 0.0;
                ny[i][j] = 0.0;
                alpha[i][j] = 1e20;
            }else{
                for (m = -1; m < 2; m++){
                    for(n = -1; n < 2; n++){
                        pg[m+1][n+1] = p[i+m][j+n];
                    }
                }
                double mx1 = 0, my1 = 0;
                double G[6], D[6][2];
                
                switch (caseno)
                {
                    case 0: //Finite Difference Method
                        mx1 = (p[i+1][j-1]-p[i-1][j-1]) +
                                2.0*(p[i+1][j]-p[i-1][j]) +
                                (p[i+1][j+1]-p[i-1][j+1]);
                        my1 = (p[i-1][j+1]-p[i-1][j-1]) +
                                2.0*(p[i][j+1]-p[i][j-1]) +
                                (p[i+1][j+1]-p[i+1][j-1]);
                        //mx1 = mx[i][j], my1 = my[i][j];
                        break;
                    case 1:  //Center Of Mass code
                        for(n = -1; n <= 1; n++){
                            for(m = -1; m <= 1; m++){
                                mx1 += pg[m+1][n+1]*((double)m);
                                my1 += pg[m+1][n+1]*((double)n);
                            }
                        }
                        break;
                    case 2 ... 3:
                        //ELVIRA CODE
                        for(m = 0; m < 3; m++)
                            G[m] = pg[m][0]+pg[m][1]+pg[m][2];//vertical sum
                        for(m = 3; m < 6; m++)
                            G[m] = pg[0][m-3]+pg[1][m-3]+pg[2][m-3]; //hortizontal sum
                        for(m = 0; m < 6; m+=3){
                            D[m][m/3] = G[m+1]-G[m];
                            D[m+1][m/3] =(G[m+2]-G[m])*0.5;
                            D[m+2][m/3] = G[m+2]-G[m+1];
                        }
                        for (m = 0; m < 3; m++){
                            D[m][1] = sign(D[4][1]);
                            D[m+3][0] = sign(D[1][0]);
                        }
                        double minerror = 10; int index;  double var;
                        for(n = 0; n < 6; n++){
                            double error = findError(atan2(D[n][1],D[n][0]));
                            if(error < minerror){
                                minerror = error;
                                index = n;
                            }
                        }
                        mx1 = D[index][0];  my1 = D[index][1];
                        if (caseno == 3){ //LVIRA Code
                            var = findmin(atan2(my1,mx1));
                            //minerror = nmsimplex(&MinDat, findError, var, 1);
                            mx1 = cos(var); my1 = sin(var);
                        }
                        break;
                }
                double nrm = fabs(mx1)+fabs(my1);
                nx[i][j] = mx1/nrm;
                ny[i][j] = my1/nrm;
                alpha[i][j] = alphacomp(1-p[i][j],nx[i][j],ny[i][j]);
            }
            if ((p[i][j]<1e-8) || (p[i][j]>(1-1e-8))){
                nx[i][j] = 0.0;
                ny[i][j] = 0.0;
                alpha[i][j] = 1e20;
            }
            alpha[i][j] = alpha[i][j] + min(0.0,nx[i][j]) + min(0.0,ny[i][j]); 
            double uthres = 1.0-1.0e-8;
            if (!BC){
                int unif = (p[i+1][j+1]>uthres) && (p[i][j+1]>uthres) &&
                    (p[i-1][j+1]>uthres) && (p[i+1][j]>uthres) &&
                    (p[i][j]>uthres) && (p[i-1][j]>uthres) &&
                    (p[i+1][j-1]>uthres) &&
                    (p[i][j-1]>uthres) && (p[i-1][j-1]>uthres);
                if (unif) alpha[i][j] = 2e19;
            }
        }

    }

    
    //free(mx); free(my); 
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Macros for the ouput and input arguments */
    double *B,  *C, *D, *psin, dx, dy, dt;
    int M, N, i, j, k, num;
    M = mxGetScalar(prhs[0]); /* Get the dimensions of A */
    N = mxGetScalar(prhs[1]);
    dx =  mxGetScalar(prhs[2]);
    dy =  dx;

    psin = mxGetPr(prhs[3]);
    caseno = mxGetScalar(prhs[4]);caseno--;
    
    plhs[0] = mxCreateDoubleMatrix(M+4, N+4, mxREAL); /* Create the output matrix */
    plhs[1] = mxCreateDoubleMatrix(M+4, N+4, mxREAL); /* Create the output matrix */
    plhs[2] = mxCreateDoubleMatrix(M+4, N+4, mxREAL); /* Create the output matrix */

    B = mxGetPr(plhs[0]);
    C =  mxGetPr(plhs[1]);
    D =  mxGetPr(plhs[2]);

    double (*psi)[N+4]; psi= malloc((M+4) * sizeof *psi );
    double (*nx )[N+4]; nx = malloc((M+4) * sizeof *nx );
    double (*ny )[N+4]; ny = malloc((M+4) * sizeof *ny );
    double (*alpha )[N+4]; alpha = malloc((M+4) * sizeof *alpha );

    
    for(j = 0; j < N+4; j++){
        for(i = 0; i < M+4; i++){
            psi[i][j] = psin[i + (M+4)*j];
        }
    }
    /***********************************************************
     *Begin Computation
     ***********************************************************/
    normalvec(M, N, dx, dy, psi, nx, ny, alpha);
    /*nx and ny are not normal but sum vectors*/
    
    /***********************************************************
     *End Computation
     ***********************************************************/
    for(i = 0; i < M+4; i++){
        for(j = 0; j < N+4; j++){
            B[i + (M+4)*j] = nx[i][j];
            C[i + (M+4)*j] = ny[i][j];
            D[i + (M+4)*j] = alpha[i][j];
        }
    }
    free(psi); free(nx); free(ny); free(alpha);
    return;
}

double eps(double x)
{
  double r, absxk;
  int exponent;
  absxk = fabs(x);
    if (absxk <= 2.2250738585072014E-308) {
      r = 4.94065645841247E-324;
    } else {
      frexp(absxk, &exponent);
      r = ldexp(1.0, exponent - 53);
    }
  return r;
}



double findmin(double x0)
{
    double x, fv[2],v[2], xe,xr,d0,fxr, d1, fxe, fvt[2];
    int lastCol,itercount,fun_evals,firstCol;
    signed char idx[2], idxb[2];
    bool doShrink,exitg1,p,guard1 = false;
    for (lastCol = 0; lastCol < 2; lastCol++) {
        v[lastCol] = x0;
    }
    
    fv[0] = findError(x0);
    if (x0 != 0.0) {
        v[1] = 1.05 * x0;
    } else {
        v[1] = 0.00025;
    }
    
    fv[1] = findError(v[1]);
    if ((fv[0] <= fv[1])) {
        idx[0] = 1;
        idx[1] = 2;
    } else {
        idx[0] = 2;
        idx[1] = 1;
    }
    
    doShrink = false;
    itercount = 1;
    fun_evals = 2;
    lastCol = idx[1] - 1;
    firstCol = idx[0] - 1;
    exitg1 = false;
    while ((!exitg1) && ((fun_evals < 200) && (itercount < 200))) {
        xe = 0.0;
        xr = fabs(fv[idx[0] - 1] - fv[idx[1] - 1]);
        if (xr > 0.0) {
            xe = xr;
        }
        
        xr = 10.0 * eps(fv[idx[0] - 1]);
        xr = 10.0 * eps(xr);
        if ((0.00001 > xr)) {
            d0 = 0.00001;
        } else {
            d0 = xr;
        }
        
        if (xe > d0) {
            p = false;
        } else {
            xe = 0.0;
            xr = fabs(v[idx[0] - 1] - v[idx[1] - 1]);
            if (xr > 0.0) {
                xe = xr;
            }
            
            xr = 10.0 * eps(v[idx[0] - 1]);
            if ((tol > xr)) {
                d1 = tol;
            } else {
                d1 = xr;
            }
            
            p = (xe <= d1);
        }
        
        if (!p) {
            xr = 2.0 * v[firstCol] - v[lastCol];
            fxr = findError(xr);
            fun_evals++;
            if (fxr < fv[idx[0] - 1]) {
                xe = 3.0 * v[firstCol] - 2.0 * v[lastCol];
                fxe = findError(xe);
                fun_evals++;
                if (fxe < fxr) {
                    v[lastCol] = xe;
                    fv[idx[1] - 1] = fxe;
                } else {
                    v[lastCol] = xr;
                    fv[idx[1] - 1] = fxr;
                }
            } else if (fxr < fv[idx[0] - 1]) {
                v[lastCol] = xr;
                fv[idx[1] - 1] = fxr;
            } else {
                guard1 = false;
                if (fxr < fv[idx[1] - 1]) {
                    x = 1.5 * v[firstCol] - 0.5 * v[lastCol];
                    xr = findError(x);
                    fun_evals++;
                    if (xr <= fxr) {
                        v[lastCol] = x;
                        fv[idx[1] - 1] = xr;
                    } else {
                        doShrink = true;
                        guard1 = true;
                    }
                } else {
                    x = 0.5 * v[firstCol] + 0.5 * v[lastCol];
                    xr = findError(x);
                    fun_evals++;
                    if (xr < fv[idx[1] - 1]) {
                        v[lastCol] = x;
                        fv[idx[1] - 1] = xr;
                    } else {
                        doShrink = true;
                        guard1 = true;
                    }
                }
                
                if (guard1) {
                    v[idx[1] - 1] = v[firstCol] + 0.5 * (v[idx[1] - 1] - v[firstCol]);
                    fv[idx[1] - 1] = findError(v[idx[1] - 1]);
                    fun_evals++;
                }
            }
            
            if (doShrink) {
                for (lastCol = 0; lastCol < 2; lastCol++) {
                    fvt[lastCol] = fv[idx[lastCol] - 1];
                    idxb[lastCol] = idx[lastCol];
                }
                
                if ((fvt[0] <= fvt[1])) {
                    idx[0] = 1;
                    idx[1] = 2;
                } else {
                    idx[0] = 2;
                    idx[1] = 1;
                }
                
                for (lastCol = 0; lastCol < 2; lastCol++) {
                    idx[lastCol] = idxb[idx[lastCol] - 1];
                }
                
                doShrink = false;
            } else {
                if (fv[idx[1] - 1] < fv[idx[0] - 1]) {
                    lastCol = idx[1];
                    idx[1] = idx[0];
                    idx[0] = (signed char)lastCol;
                }
            }
            
            itercount++;
            lastCol = idx[1] - 1;
            firstCol = idx[0] - 1;
        } else {
            exitg1 = true;
        }
    }
    //printf("%i %i\n",itercount,fun_evals);
    return v[idx[0] - 1];
}