#include <mex.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
const int INSIDE = 0; // 0000
const int LEFT = 1;   // 0001
const int RIGHT = 2;  // 0010
const int BOTTOM = 4; // 0100
const int TOP = 8;    // 1000
int ComputeOutCode(double x, double y, double xmin, double xmax, double ymin, double ymax){
    int code;
    code = INSIDE;          // initialised as being inside of [[clip window]]
    if (x < xmin)           // to the left of clip window
        code |= LEFT;
    else if (x > xmax)      // to the right of clip window
        code |= RIGHT;
    if (y < ymin)           // below the clip window
        code |= BOTTOM;
    else if (y > ymax)      // above the clip window
        code |= TOP;
    return code;
}
int CohenSutherlandLineClip(double *x0, double *y0, double *x1, double *y1, double xmin, double xmax, double ymin, double ymax){
    
    int outcode0 = ComputeOutCode(*x0, *y0, xmin, xmax, ymin, ymax);
    int outcode1 = ComputeOutCode(*x1, *y1, xmin, xmax, ymin, ymax);
    int accept = 0;
    
    while (1) {
        if (!(outcode0 | outcode1)) {
            accept = 1;
            break;
        } else if (outcode0 & outcode1) {
            break;
        } else {
            double x, y;
            
            int outcodeOut = outcode0 ? outcode0 : outcode1;
            
            if (outcodeOut & TOP) {           // point is above the clip window
                x = *x0 + (*x1 - *x0) * (ymax - *y0) / (*y1 - *y0);
                y = ymax;
            } else if (outcodeOut & BOTTOM) { // point is below the clip window
                x = *x0 + (*x1 - *x0) * (ymin - *y0) / (*y1 - *y0);
                y = ymin;
            } else if (outcodeOut & RIGHT) {  // point is to the right of clip window
                y = *y0 + (*y1 - *y0) * (xmax - *x0) / (*x1 - *x0);
                x = xmax;
            } else if (outcodeOut & LEFT) {   // point is to the left of clip window
                y = *y0 + (*y1 - *y0) * (xmin - *x0) / (*x1 - *x0);
                x = xmin;
            }
            
            if (outcodeOut == outcode0) {
                *x0 = x;
                *y0 = y;
                outcode0 = ComputeOutCode(*x0, *y0, xmin, xmax, ymin, ymax);
            } else {
                *x1 = x;
                *y1 = y;
                outcode1 = ComputeOutCode(*x1, *y1, xmin, xmax, ymin, ymax);
            }
        }
    }
    return accept;
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    srand(554);
    /* Macros for the ouput and input arguments */
    double *B, *U, *V, *mxin,  *myin, *a1in, dx, dy, *Xpts, *Ypts;
    int M, N, i, j, e = 0, ind;
    M = mxGetScalar(prhs[0]); /* Get the dimensions of A */
    N = mxGetScalar(prhs[1]);
    
    mxin = mxGetPr(prhs[2]);
    myin = mxGetPr(prhs[3]);
    a1in = mxGetPr(prhs[4]);

    dx = mxGetScalar(prhs[5]); 
    
    plhs[0] = mxCreateDoubleMatrix(2*M*N,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(2*M*N,1, mxREAL);
    
    
    Xpts = mxGetPr(plhs[0]);
    Ypts = mxGetPr(plhs[1]);
    
    for(i = 2; i < M+2; i++){//for(i = M + 1; i > 1; i--)
        for(j = 2; j < N + 2; j++){
            ind = i+(M+4)*(j);
            double alpha = a1in[ind];
            if (alpha < 1.0e15){
                double mx = mxin[ind];
                double my = myin[ind];
                double box[4][2] = {{-.5,-.5},{.5,-.5},{.5,.5},{-.5,.5}};
                double nrm = sqrt(mx*mx+my*my);
                double m1 = mx/nrm, m2 = my/nrm, a = (alpha-mx/2-my/2)/nrm;
                double x0 = a*m1+m2, y0 = a*m2-m1, x1 = a*m1-m2, y1 = a*m2+m1;
                int accept = CohenSutherlandLineClip(&x0, &y0, &x1, &y1, -.5, .5, -.5, .5);
                double x15 = ((double)i-1.5)*dx; double y15 = ((double)j-1.5)*dx;
                Xpts[e] = x15+dx*x0; Ypts[e] = y15+dx*y0; e++;
                Xpts[e] = x15+dx*x1; Ypts[e] = y15+dx*y1; e++;
                Xpts[e] = mxGetNaN();Ypts[e] = mxGetNaN(); e++;
            }
        }
    }
    while (e<2*M*N){
        Xpts[e] = mxGetNaN() ;Ypts[e] = mxGetNaN();
        e++;
    }
    return;
}