#include <mex.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include "shape.h"

const int INSIDE = 0; // 0000
const int LEFT = 1;   // 0001
const int RIGHT = 2;  // 0010
const int BOTTOM = 4; // 0100
const int TOP = 8;    // 1000


void solve(double* a, double* b, double* c, double* d, int s, int n) {
    n--;int i;
    c[s] /= b[s];
    d[s] /= b[s];
    for (i = s+1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }
    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);
    for (i = n; i-- > s;) d[i] -= c[i]*d[i+1];
}

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

void subComputePolygon(double mx, double my, double alpha, poly p1, double sx, double sy){
  double box[4][2] = {{-.5,-.5},{.5,-.5},{.5,.5},{-.5,.5}};
  int i = 0; vec_t tmp;
  if (alpha<1e19){
    double nrm = sqrt(mx*mx+my*my);
    double m1 = mx/nrm, m2 = my/nrm, a = (alpha-mx/2-my/2)/nrm;
    double x0 = a*m1+m2, y0 = a*m2-m1, x1 = a*m1-m2, y1 = a*m2+m1;
    //double x0 = a*m1+m2/2, y0 = a*m2-m1/2, x1 = a*m1-m2/2, y1 = a*m2+m1/2;
    int ret = 0, i1;
    if (isnan(x0)) printf("exists1\n");
    int accept = CohenSutherlandLineClip(&x0, &y0, &x1, &y1, -.5, .5, -.5, .5);
    if (isnan(x0)) printf("exists\n");
    
    
    
    //first points
    if (accept){
        vnew(x0+sx, y0+sy, &tmp);
        poly_append(p1, &tmp);
    }
    //all points
    if (accept) ret = ceil(atan2 (y0,x0) * 2/3.1415926535+1.5);
    for (i = ret; i<4+ret; i++){
        i1 = i%4;
        if (((m1*box[i1][0]+m2*box[i1][1])>a) || alpha>1.0e19){
            vnew(box[i1][0]+sx,box[i1][1]+sy,&tmp);
            poly_append(p1, &tmp);
        }
    }
    //last points
    if (accept){
        vnew(x1+sx, y1+sy, &tmp);
        poly_append(p1, &tmp);
    }
  }
  else{
    for (i = 0; i<4; i++){
            vnew(box[i][0]+sx,box[i][1]+sy,&tmp);
            poly_append(p1, &tmp);
    }
  }
}

void computePolygon(double mx, double my, double alpha, double psi, poly p1, double* X, double* Y){
    double nrm = sqrt(mx*mx+my*my);
    double m1 = mx/nrm, m2 = my/nrm, a = (alpha-mx/2-my/2)/nrm;
    double x0 = a*m1+m2, y0 = a*m2-m1, x1 = a*m1-m2, y1 = a*m2+m1;
    //double x0 = a*m1+m2/2, y0 = a*m2-m1/2, x1 = a*m1-m2/2, y1 = a*m2+m1/2;
    int ret = 0, i1;
    if (isnan(x0)) printf("exists1\n");
    int accept = CohenSutherlandLineClip(&x0, &y0, &x1, &y1, -.5, .5, -.5, .5);
    if (isnan(x0)) printf("exists\n");
    double box[8][2] = {{-.5,-.5},{0,-.5},{.5,-.5},{.5,0},{.5,.5},{0,.5},{-.5,.5},{-.5,0}};
    int i = 0; vec_t tmp;
    double N[8], xn0, yn0, xn1, yn1, xr0, yr0, xr1, yr1;
    //first points
    if (accept){
        xn0 = 0; yn0 = 0;
        shape(x0, y0, N);
        for (i = 0; i<8; i++){
            xn0+=X[i]*N[i]; yn0+=Y[i]*N[i];
        }
        vnew(xn0, yn0, &tmp);
        poly_append(p1, &tmp);
    }
    //all points
    if (accept) ret = ceil(atan2 (y0,x0) * 4/3.1415926535+3);
    for (i = ret; i<8+ret; i++){
        i1 = i%8;
        if (((m1*box[i1][0]+m2*box[i1][1])>a) || alpha>1.0e19){
            vnew(X[i1],Y[i1],&tmp);
            poly_append(p1, &tmp);
        }
    }
    //last points
    if (accept){
        xn1 = 0; yn1 = 0;
        shape(x1, y1, N);
        for (i = 0; i<8; i++){
            xn1+=X[i]*N[i]; yn1+=Y[i]*N[i];
        }
        vnew(xn1, yn1, &tmp);
        poly_append(p1, &tmp);
    }
    if (accept){
        //area correction point
        double ncase = 4.0*mx*my;
        int samesign = ((x0>0)&(y0>0)&(x1>0)&(y1>0)) || ((x0<0)&(y0<0)&(x1<0)&(y1<0));
        int diffsign = ((x0>0)&(y0<0)&(x1>0)&(y1<0)) || ((x0<0)&(y0>0)&(x1<0)&(y1>0));
        
        //printf("\n\n%i %i %f\n\n",samesign,diffsign,ncase);
        //(x0 y0 x1 y1>0) or (x0 y0 x1 y1<0) --> 0--4
        //(x0 x1<0 y0 y1>0) or (x0 x1>0 y0 y1<0) --> 2--6
        if ((ncase>0.9) || samesign){
            xr0 = X[0]; yr0 = Y[0]; xr1 = X[4]; yr1 = Y[4];
        }
        else if ((ncase<-0.9) || diffsign){
            xr0 = X[2]; yr0 = Y[2]; xr1 = X[6]; yr1 = Y[6];
        }
        else{
            if (fabs(mx)>fabs(my)){
                xr0 = X[7]; yr0 = Y[7]; xr1 = X[3]; yr1 = Y[3];
            }
            else{
                xr0 = X[1]; yr0 = Y[1]; xr1 = X[5]; yr1 = Y[5];
            }
        }
        //solve point alloc
        double An = psi-polygonArea(p1);
        //printf("%f\n",An);
        double ad[] = {0,yr0-yr1}; double bd[]={yn0-yn1,xr1-xr0};
        double cd[] = {xn1-xn0}; double dv[] = {xn1*yn0-xn0*yn1+2.0*An,
        -(xr0-xr1)*yr0-(yr1-yr0)*xr0}; solve(ad,bd,cd,dv,0,2);
        //printf("%f %f\n",dv[0],dv[1]);
        vnew(dv[0],dv[1],&tmp);
        poly_append(p1, &tmp);
    }
}

double polygonArea(poly p1){
    double area = 0.0;
    int i, j = p1->len-1;
    for (i=0; i<p1->len; i++) {
        area -= (p1->v[j].x+p1->v[i].x)*(p1->v[j].y-p1->v[i].y); j = i;
    }
    return area*.5;
}

void shape(double p1, double p2, double *N){
    N[0] =  (fabs(p1)-p1)*(fabs(p2)-p2);
    N[1] =  (-2.0*fabs(p1)+1.0)*(fabs(p2)-p2);
    N[2] =  (fabs(p1)+p1)*(fabs(p2)-p2);
    N[3] =  (fabs(p1)+p1)*(-2.0*fabs(p2)+1.0);
    N[4] =  (fabs(p1)+p1)*(fabs(p2)+p2);
    N[5] =  (-2.0*fabs(p1)+1.0)*(fabs(p2)+p2);
    N[6] =  (fabs(p1)-p1)*(fabs(p2)+p2);
    N[7] =  (fabs(p1)-p1)*(-2.0*fabs(p2)+1.0);
}

double dot(vec a, vec b) {
    return a->x * b->x + a->y * b->y;
}

double cross(vec a, vec b) {
	return a->x * b->y - a->y * b->x;
}

vec vsub(vec a, vec b, vec res) {
	res->x = a->x - b->x;
	res->y = a->y - b->y;
	return res;
}

vec vnew(double a, double b, vec res){
	res->x = a;
	res->y = b;
	return res;
}

/* tells if vec c lies on the left side of directed edge a->b
 * 1 if left, -1 if right, 0 if colinear
 */
int left_of(vec a, vec b, vec c) {
	vec_t tmp1, tmp2;
	double x;
	vsub(b, a, &tmp1);
	vsub(c, b, &tmp2);
	x = cross(&tmp1, &tmp2);
	return x < 0 ? -1 : x > 0;
}

int line_sect(vec x0, vec x1, vec y0, vec y1, vec res) {
	vec_t dx, dy, d;
	vsub(x1, x0, &dx);
	vsub(y1, y0, &dy);
	vsub(x0, y0, &d);
	/* x0 + a dx = y0 + b dy ->
	   x0 X dx = y0 X dx + b dy X dx ->
	   b = (x0 - y0) X dx / (dy X dx) */
	double dyx = cross(&dy, &dx);
	if (!dyx)
		return 0;
	dyx = cross(&d, &dx) / dyx;
	if (dyx <= 0 || dyx >= 1)
		return 0;

	res->x = y0->x + dyx * dy.x;
	res->y = y0->y + dyx * dy.y;
	return 1;
}

/* === polygon stuff === */

poly poly_new() {
	return (poly)calloc(1, sizeof(poly_t));
}

void poly_free(poly p) {
	free(p->v);
	free(p);
}

void poly_append(poly p, vec v) {
	if (p->len >= p->alloc) {
		p->alloc *= 2;
		if (!p->alloc)
			p->alloc = 4;
		p->v = (vec)realloc(p->v, sizeof(vec_t) * p->alloc);
	}
	p->v[p->len++] = *v;
}

/* this works only if all of the following are true:
 *   1. poly has no colinear edges;
 *   2. poly has no duplicate vertices;
 *   3. poly has at least three vertices;
 *   4. poly is convex (implying 3).
*/
int poly_winding(poly p) {
	return left_of(p->v, p->v + 1, p->v + 2);
}

void poly_edge_clip(poly sub, vec x0, vec x1, int left, poly res) {
	int i, side0, side1;
	vec_t tmp;
	vec v0 = sub->v + sub->len - 1, v1;
	res->len = 0;

	side0 = left_of(x0, x1, v0);
	if (side0 != -left)
		poly_append(res, v0);

	for (i = 0; i < sub->len; i++) {
		v1 = sub->v + i;
		side1 = left_of(x0, x1, v1);
		if (side0 + side1 == 0 && side0)
			/* last point and current straddle the edge */
			if (line_sect(x0, x1, v0, v1, &tmp))
				poly_append(res, &tmp);
		if (i == sub->len - 1)
			break;
		if (side1 != -left)
			poly_append(res, v1);
		v0 = v1;
		side0 = side1;
	}
}

poly poly_clip(poly sub, poly clip) {
	int i;
	poly p1 = poly_new(), p2 = poly_new(), tmp;

	int dir = poly_winding(clip);
	poly_edge_clip(sub, clip->v + clip->len - 1, clip->v, dir, p2);
	for (i = 0; i < clip->len - 1; i++) {
		tmp = p2;
		p2 = p1;
		p1 = tmp;
		if (p1->len == 0) {
			p2->len = 0;
			break;
		}
		poly_edge_clip(p1, clip->v + i, clip->v + i + 1, dir, p2);
	}

	poly_free(p1);
	return p2;
}