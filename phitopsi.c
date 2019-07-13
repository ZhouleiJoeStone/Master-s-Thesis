#include "mex.h"
#include "math.h"
#include "stdint.h"

double min(double a, double b) {return ((a > b) ? b : a);}
double max(double a, double b) {return ((a > b) ? a : b);}
double sign(double x) {return ((x > 0) ? 1.0 : ((x < 0) ? -1.0 : 0.0));}
double step(double x) {return ((x > 0) ? 1.0 : ((x < 0) ?  0.0 : 0.5));}
int caseno;
double norm(double x, double y){
  return sqrt(x*x+y*y);
}

double zalezak(double x, double y){
  double R = sqrt(pow(x-.5,2)+pow(y-.75,2));
  double phi = .15-R;
  double bottom = .75-.15*cos(asin(.025/.15));
  double dfc = phi;
  if ((y > .85) && (R < .15))
    phi = min(dfc,y-.85);
  if ((y < .85) && (x < .475) && (R < .15))
    phi = min(dfc,.475-x);
  if ((y > .85) && (x < .475) && (R < .15))
    phi = min(dfc,sqrt(pow(.475-x,2)+pow(y-.85,2)));
  if ((y < .85) && (x > .525) && (R < .15))
    phi = min(dfc,x-.525);
  if ((y > .85) && (x > .525) && (R < .15))
    phi = min(dfc,sqrt(pow(.525-x,2)+pow(y-.85,2)));
  if ((y < .85) && (x > .525) && (R < .15))
    phi = min(dfc,x-.525);
  if ((x > .475) && (x < .525) && (y < .85) && (y > bottom))
    phi = min(0.85-y,min(.525-x,x-.475));
  if ((x > .475) && (x < .525) && (y < bottom))
    phi = min(norm(bottom-y,x-.475),norm(bottom-y,x-.525));
  double sgn = (R<.15);
  if ((x>.475) && (x<.525) && (y<.85))
    sgn = 0;
  sgn = 2.0*sgn-1.0;
  phi = fabs(phi)*sgn;
  return phi;
}

double ellipse(double px, double py, double sx, double sy, double x, double y){
  return 1.0-sqrt(pow(1.0/sx*(x-px),2)+pow(1.0/sy*(y-py),2));
}

double f(double x, double y){
    if (caseno == 2){
    double sc = 1.0;
    return max(ellipse(0.5/sc,0.75/sc,.15/sc,.15/sc,x,y),
            ellipse(0.5/sc,-1.0/sc,10.75/sc,.05/sc,x,y));}
    else
        return zalezak(x, y);		
}

double topsi(double lb, double rb, double lt, double rt){
    if ((lb>0.0) && (rb>0.0) && (lt>0.0) && (rt>0.0)){
        return 1.0;
    }
    else if ((lb<0.0) && (rb<0.0) && (lt<0.0) && (rt<0.0)){
        return 0.0;
    }
    else{
        double a,b;
        if ((lb*lt)<0 && (lb*rb)<0){
            a = fabs(lb)/(fabs(lb)+fabs(lt));
            b = fabs(lb)/(fabs(lb)+fabs(rb));
            if (lb>0.0) return 0.5*a*b;
            else return 1.0-0.5*a*b;
        }
        else if ((rb*lb)<0 && (rb*rt)<0){
            a = fabs(rb)/(fabs(rb)+fabs(lb));
            b = fabs(rb)/(fabs(rb)+fabs(rt));
            if (rb>0.0) return 0.5*a*b;
            else return 1.0-0.5*a*b;
        }
        else if ((rt*rb)<0 && (rt*lt)<0){
            a = fabs(rt)/(fabs(rt)+fabs(rb));
            b = fabs(rt)/(fabs(rt)+fabs(lt));
            if (rt>0.0) return 0.5*a*b;
            else return 1.0-0.5*a*b;
        }
        else if ((lt*rt)<0 && (lt*lb)<0){
            a = fabs(lt)/(fabs(lt)+fabs(rt));
            b = fabs(lt)/(fabs(lt)+fabs(lb));
            if (lt>0.0) return 0.5*a*b;
            else return 1.0-0.5*a*b;
        }
        else if ((lt*lb>0) && (rt*rb>0) && (lb*rb<0) && (lt*rt<0)){
            a = fabs(lt)/(fabs(lt)+fabs(rt));
            b = fabs(lb)/(fabs(lb)+fabs(rb));
            if (lt>0.0) return 0.5*(a+b);
            else return 1.0-0.5*(a+b);
        }
        else if (!(lt*lb>0) && !(rt*rb>0) && !(lb*rb<0) && !(lt*rt<0)){
            a = fabs(rb)/(fabs(rb)+fabs(rt));
            b = fabs(lb)/(fabs(lb)+fabs(lt));
            if (rb>0.0) return 0.5*(a+b);
            else return 1.0-0.5*(a+b);
        }
    }
}
int count = 0;
double quadtree(double x1, double x2, double y1, double y2, double err){
  count++;
  if (f(x1,y1)>0 && f(x1,y2)>0 && f(x2,y1)>0 && f(x2,y2)>0)
    return (x2-x1)*(y2-y1);
  else if (f(x1,y1)<0 && f(x1,y2)<0 && f(x2,y1)<0 && f(x2,y2)<0)
    return 0;
  else{
    if ((x2-x1)*(y2-y1)>err){
      double x12 = (x1+x2)/2.0;
      double y12 = (y1+y2)/2.0;
      return quadtree(x1,x12,y1,y12,err) +  quadtree(x12,x2,y1,y12,err) +
      quadtree(x1,x12,y12,y2,err) +  quadtree(x12,x2,y12,y2,err);}
    else
      return topsi(f(x1,y1),f(x2,y1),f(x1,y2),f(x2,y2))*(x2-x1)*(y2-y1);
      }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Macros for the ouput and input arguments */
    double *psi, *phi, dx, dy, *x, *y;
    int M, N, i, j, k, num;
    M = mxGetScalar(prhs[0]); 
    N = mxGetScalar(prhs[1]);
    x = mxGetPr(prhs[2]);
    y = mxGetPr(prhs[3]);
    dx = mxGetScalar(prhs[4]); /* Get the dimensions of A */
    dy = mxGetScalar(prhs[5]);
    caseno = mxGetScalar(prhs[6]); 

    plhs[0] = mxCreateDoubleMatrix(M+4, N+4, mxREAL); /* Create the output matrix */
    
    psi = mxGetPr(plhs[0]);
    
    for(j = 0; j < N+4; j++){
        for(i = 0; i < M+4; i++){
            psi[i + (M+4)*j] = 0;
        }
    }
    /***********************************************************
     *Begin Computation
     ***********************************************************/
    double phic, lb, rb, lt, rt;
    for(j = 1; j < N+3; j++){
        for(i = 1; i < M+3; i++){
            int ind = i + (M+4)*j;
            psi[ind] = quadtree(x[ind]-dx/2.0,
                    x[ind]+dx/2.0,
                    y[ind]-dy/2.0,
                    y[ind]+dy/2.0,
                    dx*dy*1e-7)/(dx*dy);
        }
    }

    /*nx and ny are not normal but sum vectors*/
    
    /***********************************************************
     *End Computation
     ***********************************************************/
    
    return;
}
