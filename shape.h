typedef struct { double x, y; } vec_t, *vec;

double dot(vec a, vec b);
double cross(vec a, vec b);
vec vsub(vec a, vec b, vec res);
vec vnew(double a, double b, vec res);
int left_of(vec a, vec b, vec c);
int line_sect(vec x0, vec x1, vec y0, vec y1, vec res);

/* === polygon stuff === */
typedef struct {int len, alloc;	vec v;} poly_t, *poly;
poly poly_new();
void poly_free(poly p);
void poly_append(poly p, vec v);
int poly_winding(poly p);
void poly_edge_clip(poly sub, vec x0, vec x1, int left, poly res);
poly poly_clip(poly sub, poly clip);
double polygonArea(poly p1);

/* === Other Subroutines */
void shape(double p1, double p2, double *N);
void solve(double* a, double* b, double* c, double* d, int s, int n);
int ComputeOutCode(double x, double y, double xmin, double xmax, double ymin, double ymax);
int CohenSutherlandLineClip(double *x0, double *y0, double *x1, double *y1, double xmin, double xmax, double ymin, double ymax);
void computePolygon(double mx, double my, double alpha, double psi, poly p1, double* X, double* Y);
void subComputePolygon(double mx, double my, double alpha, poly p1, double sx, double sy);
