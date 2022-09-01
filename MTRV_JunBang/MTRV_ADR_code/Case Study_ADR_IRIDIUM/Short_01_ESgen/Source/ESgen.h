double mu = 398600.442;
double pi = 3.1415926535;
double d2r = 3.1415926535 / 180;
double re = 6378.137;
double minute = 60;
double hr = 3600;

double eps_xg = 0.001,	eps_bnd = 0.01;
double eps_same_t = 20,	eps_same_J = 0.0001;
double min_gridsize = 1800;

double lb_ttr, ub_ttr;
double Q1[6], Q2[6], tmax;
double Jval[3], sol_gd[10000], x_local[30000];

void delV(double *, double *, double *, double *);
void getstate(double *, double, double *, double *);
void TLAMB(int, double, double, double, int, double *);
void XLAMB(int, double, double, double, double *);
void VLAMB(double, double, double, double, double *);

void findlocal(double*, double*, double, double*);
void findlocal2(double*, double*, double, double*);
void optim(double*, double*, double*, double *);
int usrf_(long *Status, long *n, double x[],
	long *needF, long *nF, double F[],
	long *needG, long *neG, double G[],
	char *cu, long *lencu,
	long iu[], long *leniu,
	double ru[], long *lenru);

double NORM(double*);
void NORMA(double*, double*);
double DOT(double*, double*);
void CROSS(double*, double*, double*);
double MIN(double*, int);
double MAX(double*, int);
int issame(double*, double*, double*, double*);
void xupdate(double*, double*, double*);

FILE *finp;
FILE *flocal;