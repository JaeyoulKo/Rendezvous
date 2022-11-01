#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "ESgen.h"
#include "snopt.h"
#include "f2c.h"
#include "snfilewrapper.h"

int main(){
	clock_t c1, c2;
	c1 = clock();

	int ntarget, ind, i, j;
	double *Q;
	char fname[25];

	ntarget = 100;
	tmax = 7*24*3600;
	lb_ttr = 0.5*3600;
	ub_ttr = 6*3600;

	finp = fopen("target_list.txt", "r");
	Q = (double*)calloc(6 * (ntarget + 1), sizeof(double));
	for (ind = 1; ind < ntarget + 1; ind++)
		fscanf(finp, "%lf %lf %lf %lf %lf %lf\n", &(Q[6 * ind]), &(Q[6 * ind + 1]), &(Q[6 * ind + 2]), 
		&(Q[6 * ind + 3]), &(Q[6 * ind + 4]), &(Q[6 * ind + 5]));
	fclose(finp);

	mkdir(".\\Result");
	for (i = 79; i <= ntarget; i++){
		sprintf(fname, "%s%d%s", ".\\Result\\x_local", i,".txt");
		flocal = fopen(fname, "w");
		for (j = 1; j <= ntarget; j++){
			if (i == 0)
				fprintf(flocal, "%d  %d  %d  %d\n", j, 0, 0, 0);
			else if (i != j){
				printf("%d -> %d: ", i, j);
				for (ind = 0; ind < 6; ind++){
					Q1[ind] = Q[6 * i + ind];
					Q2[ind] = Q[6 * j + ind];
				}
				findlocal(Q1, Q2, tmax, x_local);
				for (ind = 0; ind < x_local[0]; ind++)
					fprintf(flocal, "%d  %lf  %lf  %lf\n", 1000 * i + j, x_local[3 * ind + 1], x_local[3 * ind + 2], x_local[3 * ind + 3]);
			}
		}
		fclose(flocal);
	}

	c2 = clock();
	printf("Runtime: %lf (sec)\n", (double)(c2-c1)/CLOCKS_PER_SEC);
	return 0;
}

void findlocal(double *Q1, double *Q2, double tmax, double *x_local){
	//x_local = [nlocal, td, ttr, J, ...]
	double n1, n2, d_td, d_ttr;
	int ngrid, ngrid_new, ngrid_td, ngrid_ttr;
	double *grid, *grid_new, *grid_td, *grid_ttr;

	double xlb[2], xub[2];
	double xg1[2], xg2[2], xg3[2], xg4[2];
	double xopt1[5], xopt2[5], xopt3[5], xopt4[5];
	double xshort[3], xlong[3];
	int nlocal = 0, nbnd, check;

	double tmp[2], tmp2[2];
	int ind, ind2;

	n1 = sqrt(mu / pow(Q1[0], 3));
	n2 = sqrt(mu / pow(Q2[0], 3));

	tmp[0] = 2 * pi / fabs(n1 - n2);
	tmp[1] = pi / n1;
	tmp2[0] = 2 * pi / n2;
	tmp2[1] = pi*fabs(n1 - n2) / n1 / n2;
	d_td = 0.9*tmp[1];
	d_ttr = 0.9*tmp2[0];

	ngrid_td = ceil(tmax / d_td);
	ngrid_ttr = ceil((ub_ttr - lb_ttr) / d_ttr);
	grid_td = (double *)calloc(ngrid_td + 1, sizeof(double));
	grid_ttr = (double *)calloc(ngrid_ttr + 1, sizeof(double));
	grid = (double *)calloc(4 * ngrid_td * ngrid_ttr, sizeof(double));

	for (ind = 0; ind < ngrid_td + 1; ind++)
		grid_td[ind] = ind*d_td;
	if (grid_td[ngrid_td] > tmax)
		grid_td[ngrid_td] = tmax;
	for (ind = 0; ind < ngrid_ttr + 1; ind++)
		grid_ttr[ind] = lb_ttr + ind*d_ttr;
	if (grid_ttr[ngrid_ttr] > ub_ttr)
		grid_ttr[ngrid_ttr] = ub_ttr;

	ngrid = 0;
	for (ind = 0; ind < ngrid_td; ind++){
		for (ind2 = 0; ind2 < ngrid_ttr; ind2++){
			grid[4 * ngrid] = grid_td[ind];
			grid[4 * ngrid + 1] = grid_td[ind + 1];
			grid[4 * ngrid + 2] = grid_ttr[ind2];
			grid[4 * ngrid + 3] = grid_ttr[ind2 + 1];
			ngrid = ngrid + 1;
		}
	}

	while (ngrid>0){
		grid_new = (double *)calloc(16 * ngrid, sizeof(double));
		ngrid_new = 0;

		for (ind = 0; ind < ngrid; ind++){
			if (grid[4 * ind] + grid[4 * ind + 2] < tmax){
				xshort[2] = 1000, xlong[2] = 1000;
				check = 2;

				tmp[0] = grid[4 * ind + 1];		tmp[1] = tmax - grid[4 * ind + 2];
				tmp2[0] = grid[4 * ind + 3];	tmp2[1] = tmax - grid[4 * ind];
				xlb[0] = grid[4 * ind];			xlb[1] = grid[4 * ind + 2];
				xub[0] = MIN(tmp, 2);			xub[1] = MIN(tmp2, 2);

				xg1[0] = xlb[0] + eps_xg;			xg1[1] = xlb[1] + eps_xg;
				xg2[0] = xub[0] - eps_xg;			xg2[1] = xub[1] - eps_xg;
				xg3[0] = xlb[0] + eps_xg;			xg3[1] = xub[1] - eps_xg;
				xg4[0] = xub[0] - eps_xg;			xg4[1] = xlb[1] + eps_xg;

				optim(xg1, xlb, xub, xopt1);
				optim(xg2, xlb, xub, xopt2);
				optim(xg3, xlb, xub, xopt3);
				optim(xg4, xlb, xub, xopt4);
				nbnd = xopt1[3] + xopt2[3] + xopt3[3] + xopt4[3];

				if (nbnd == 4)
					check = 0;
				else if (nbnd == 0 || xub[0] - xlb[0] < min_gridsize || xub[1] - xlb[1] < min_gridsize){
					check = 1;
					xupdate(xopt1, xshort, xlong);
					xupdate(xopt2, xshort, xlong);
					xupdate(xopt3, xshort, xlong);
					xupdate(xopt4, xshort, xlong);
				}

				switch (check){
				case 0:
					break;
				case 1:
					if (xshort[2] < 1000){
						x_local[3 * nlocal + 1] = xshort[0];
						x_local[3 * nlocal + 2] = xshort[1];
						x_local[3 * nlocal + 3] = xshort[2];
						nlocal = nlocal + 1;
					}
					if (xlong[2] < 1000){
						x_local[3 * nlocal + 1] = xlong[0];
						x_local[3 * nlocal + 2] = xlong[1];
						x_local[3 * nlocal + 3] = xlong[2];
						nlocal = nlocal + 1;
					}
					break;
				case 2:
					tmp[0] = 0.5*(grid[4 * ind] + grid[4 * ind + 1]);
					tmp[1] = 0.5*(grid[4 * ind + 2] + grid[4 * ind + 3]);
					grid_new[4 * ngrid_new] = grid[4 * ind];
					grid_new[4 * ngrid_new + 1] = tmp[0];
					grid_new[4 * ngrid_new + 2] = grid[4 * ind + 2];
					grid_new[4 * ngrid_new + 3] = tmp[1];

					grid_new[4 * ngrid_new + 4] = grid[4 * ind];
					grid_new[4 * ngrid_new + 5] = tmp[0];
					grid_new[4 * ngrid_new + 6] = tmp[1];
					grid_new[4 * ngrid_new + 7] = grid[4 * ind + 3];

					grid_new[4 * ngrid_new + 8] = tmp[0];
					grid_new[4 * ngrid_new + 9] = grid[4 * ind + 1];
					grid_new[4 * ngrid_new + 10] = grid[4 * ind + 2];
					grid_new[4 * ngrid_new + 11] = tmp[1];
					
					grid_new[4 * ngrid_new + 12] = tmp[0];
					grid_new[4 * ngrid_new + 13] = grid[4 * ind + 1];
					grid_new[4 * ngrid_new + 14] = tmp[1];
					grid_new[4 * ngrid_new + 15] = grid[4 * ind + 3];
					ngrid_new = ngrid_new + 4;
					break;
				}
			}
		}
		free(grid);
		ngrid = ngrid_new;
		if (ngrid == 0)
			break;
		else{
			grid = (double *)calloc(4 * ngrid, sizeof(double));
			for (ind = 0; ind < 4 * ngrid; ind++){
				grid[ind] = grid_new[ind];
			}
			free(grid_new);
		}
	}
	printf("%d local solutions were found\n", nlocal);
	x_local[0] = nlocal;
	return;
}

void optim(double *x0, double *x_lb, double *x_ub, double *opt_result){
	//opt_result: [td, ttr, J, isbnd, thetaf<pi]
	int n = 2, nF = 2;

	int isbnd = 0;
	int iPrint = -1, iSumm = -1, Start = 0;
	int lenA = n*nF, lenG = n*nF;
	int lencw = 1000, leniw = 1000, lenrw = 1000;
	int nxname = 1, nFname = 1, npname = 1;
	int INFO, neA, neG, mincw, miniw, minrw, nInf, nS, strOpt_len;
	int ObjRow = 1, ObjAdd = 0, DerOpt = 0;
	double sInf;
	char strOpt[50], Prob[10];

	char *cw, *xnames, *Fnames;
	int *iw, *iAfun, *jAvar, *iGfun, *jGvar, *xstate, *Fstate;
	double *rw, *A, *x, *xlow, *xupp, *xmul, *F, *Flow, *Fupp, *Fmul;

	rw = (double *)calloc(lenrw, sizeof(double));
	iw = (int *)calloc(leniw, sizeof(int));
	cw = (char *)calloc(lencw * 8, sizeof(char));

	iAfun = (int *)calloc(lenA, sizeof(int));
	jAvar = (int *)calloc(lenA, sizeof(int));
	iGfun = (int *)calloc(lenG, sizeof(int));
	jGvar = (int *)calloc(lenG, sizeof(int));

	xstate = (int *)calloc(n, sizeof(int));
	Fstate = (int *)calloc(nF, sizeof(int));
	A = (double *)calloc(lenA, sizeof(double));

	x = (double *)calloc(n, sizeof(double));
	xlow = (double *)calloc(n, sizeof(double));
	xupp = (double *)calloc(n, sizeof(double));
	xmul = (double *)calloc(n, sizeof(double));

	F = (double *)calloc(nF, sizeof(double));
	Flow = (double *)calloc(nF, sizeof(double));
	Fupp = (double *)calloc(nF, sizeof(double));
	Fmul = (double *)calloc(nF, sizeof(double));

	xnames = (char *)calloc(nxname * 8, sizeof(char));
	Fnames = (char *)calloc(nFname * 8, sizeof(char));

	for (int ind = 0; ind < n; ind++){
		xlow[ind] = x_lb[ind];
		xupp[ind] = x_ub[ind];
		x[ind] = x0[ind];
	}
	Flow[0] = 0;  Fupp[0] = 100;
	Flow[1] = 0;  Fupp[1] = tmax;

	sprintf(strOpt, "%s", "Derivative option = 0");
	strOpt_len = strlen(strOpt);

	sninit_(&iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw, 8 * lencw);

	snjac_(&INFO, &nF, &n, usrf_, iAfun, jAvar, &lenA, &neA, A, iGfun, jGvar, &lenG, &neG,
		x, xlow, xupp, &mincw, &miniw, &minrw, cw, &lencw, iw, &leniw, rw, &lenrw,
		cw, &lencw, iw, &leniw, rw, &lenrw,	8 * lencw, 8 * lencw);

	snset_(strOpt, &iPrint, &iSumm, &INFO,
		cw, &lencw, iw, &leniw, rw, &lenrw, strOpt_len, 8 * lencw);

	snopta_(&Start, &nF, &n, &nxname, &nFname, &ObjAdd, &ObjRow, Prob, usrf_, iAfun, jAvar,
		&lenA, &neA, A, iGfun, jGvar, &lenG, &neG, xlow, xupp, xnames, Flow, Fupp, Fnames,
		x, xstate, xmul, F, Fstate, Fmul, &INFO, &mincw, &miniw, &minrw, &nS, &nInf, &sInf,
		cw, &lencw, iw, &leniw, rw, &lenrw, cw, &lencw, iw, &leniw, rw, &lenrw,
		8 * npname, 8 * nxname, 8 * nFname, 8 * lencw, 8 * lencw); 

	if (fabs(x[0] - xlow[0]) < eps_bnd && xlow[0] != 0)
		isbnd = 1;
	else if (fabs(x[1] - xlow[1]) < eps_bnd && xlow[1] != lb_ttr)
		isbnd = 1;
	else if (fabs(x[0] - xupp[0]) < eps_bnd)
		isbnd = 1;
	else if (fabs(x[1] - xupp[1]) < eps_bnd && xupp[1] != ub_ttr)
		isbnd = 1;
	else if (x[0] < 0 || (x[0] < eps_bnd && fabs(F[1] - tmax) < eps_bnd))
		isbnd = 1;

	delV(Q1, Q2, x, Jval);
	opt_result[0] = x[0];
	opt_result[1] = x[1];
	opt_result[2] = Jval[0];
	opt_result[3] = isbnd;
	opt_result[4] = Jval[1];

	free(rw);	free(iw);	free(cw);
	return;
}

int usrf_(long *Status, long *n, double x[], long *needF, long *nF, double F[], long *needG, long *neG, 
		  double G[], char *cu, long *lencu, long iu[], long *leniu, double ru[], long *lenru)
{
	int ind0;
	delV(Q1, Q2, x, Jval);
	F[0] = Jval[0];
	F[1] = x[0] + x[1];

	return 0;
}

void delV(double *Q1, double *Q2, double *x, double *Jval){
	//Jval = [J, thetaf<pi]
	if (x[1] < eps_xg){
		Jval[0] = 1e22;
		Jval[1] = 1;
		return;
	}
	double r1[3], v1[3], r2[3], v2[3];
	double tf, rm1, rm2, thetaf;
	double er1[3], er2[3], et1[3], et2[3];
	double vr1, vt1, vr2, vt2, delv1[3], delv2[3];
	double *J1, *J2, J1min, J2min;
	int ind1, ind2;
	double tmp, tmp1[3], tmp2[3];

	getstate(Q1, x[0], r1, v1);
	getstate(Q2, (x[0]+x[1]), r2, v2);
	tf = x[1];

	rm1 = NORM(r1);
	rm2 = NORM(r2);
	tmp = DOT(r1, r2) / rm1 / rm2;
	if (tmp > 1)
		tmp = 1;
	else if (tmp < -1)
		tmp = -1;
	thetaf = acos(tmp);

	NORMA(r1, er1);
	NORMA(r2, er2);
	CROSS(r1, r2, tmp1);
	CROSS(tmp1, r1, tmp2);
	NORMA(tmp2, et1);
	CROSS(tmp1, r2, tmp2);
	NORMA(tmp2, et2);
	VLAMB(rm1, rm2, thetaf / d2r, tf, sol_gd);
	ind1 = sol_gd[0];
	J1 = (double *)calloc(ind1, sizeof(double));
	for (int i = 0; i < ind1; i++){
		vr1 = sol_gd[5 * i + 2];
		vt1 = sol_gd[5 * i + 3];
		vr2 = sol_gd[5 * i + 4];
		vt2 = sol_gd[5 * i + 5];
		for (int j = 0; j < 3; j++){
			delv1[j] = vr1 * er1[j] + vt1 * et1[j] - v1[j];
			delv2[j] = v2[j] - vr2 * er2[j] - vt2 * et2[j];
		}
		J1[i] = NORM(delv1) + NORM(delv2);
	}
	J1min = MIN(J1, ind1);

	VLAMB(rm1, rm2, 360 - thetaf / d2r, tf, sol_gd);
	ind2 = sol_gd[0];
	J2 = (double *)calloc(ind2, sizeof(double));
	for (int i = 0; i < ind2; i++){
		vr1 = sol_gd[5 * i + 2];
		vt1 = sol_gd[5 * i + 3];
		vr2 = sol_gd[5 * i + 4];
		vt2 = sol_gd[5 * i + 5];
		for (int j = 0; j < 3; j++){
			delv1[j] = vr1 * er1[j] - vt1 * et1[j] - v1[j];
			delv2[j] = v2[j] - vr2 * er2[j] + vt2 * et2[j];
		}
		J2[i] = NORM(delv1) + NORM(delv2);
	}
	J2min = MIN(J2, ind2);
	free(J1);
	free(J2);

	if (J1min < J2min){
		Jval[0] = J1min;
		Jval[1] = 1;
	}
	else{
		Jval[0] = J2min;
		Jval[1] = -1;
	}
}

void getstate(double *Q, double t1, double *X, double *Xdot){
	double a, e, i, ome, w, nu;
	double ee, nu0, M0, M, E, E_new;
	double p, r, h;
	double err, tmp1, tmp2;

	a = Q[0];
	e = Q[1];
	i = Q[2] * d2r;
	ome = Q[3] * d2r;
	w = Q[4] * d2r;
	nu0 = Q[5] * d2r;

	ee = e*e;
	//M0 = nu0 - 2 * e*sin(nu0) + (3 * ee / 4 + ee*ee / 8)*sin(2 * nu0) - e*ee*sin(3 * nu0) / 3;
	M0 = nu0;
	M = M0 + sqrt(mu / a / a / a)*t1;
	if (M < 0)
		M = M + 2 * pi;
	else if (M > 2 * pi)
		M = fmod(M, 2 * pi);

	if (e == 0)
		nu = M;
	else{
		if (M > pi)
			E = M - e;
		else
			E = M + e;
		while (1){
			E_new = E + (M - E + e*sin(E)) / (1 - e*cos(E));
			err = fabs(E_new - E);
			if (err < 1e-8)
				break;
			else
				E = E_new;
		}
		tmp1 = sin(E)*sqrt(1 - ee) / (1 - e*cos(E));
		tmp2 = (cos(E) - e) / (1 - e*cos(E));
		nu = atan2(tmp1, tmp2);
	}
	p = a*(1 - ee);
	r = p / (1 + e*cos(nu));
	h = sqrt(mu*p);

	tmp1 = cos(w + nu);
	tmp2 = sin(w + nu);
	X[0] = r*(cos(ome)*tmp1 - sin(ome)*tmp2*cos(i));
	X[1] = r*(sin(ome)*tmp1 + cos(ome)*tmp2*cos(i));
	X[2] = r*tmp2*sin(i);
	Xdot[0] = X[0] * h*e*sin(nu) / r / p - h*(cos(ome)*tmp2 + sin(ome)*tmp1*cos(i)) / r;
	Xdot[1] = X[1] * h*e*sin(nu) / r / p - h*(sin(ome)*tmp2 - cos(ome)*tmp1*cos(i)) / r;
	Xdot[2] = X[2] * h*e*sin(nu) / r / p + h*tmp1*sin(i) / r;
}

void TLAMB(int m, double q, double qsqfm1, double x, int N, double *tval){
	//tval = [t,dt,ddt,dddt]
	double t, dt, ddt, dddt;
	int LM1, L1, L2, L3;
	double u, y, z, qq, xx, qx, f, g, fg1, t_old;
	double A, B, AA, BB;
	double tmp1, tmp2;

	t = -100;
	LM1 = N == -1;	L1 = N >= 1;	L2 = N >= 2;	L3 = N == 3;
	qq = q*q;
	xx = x*x;

	u = 1 - xx;
	if (!LM1){ dt = 0;	ddt = 0; dddt = 0; }
	if (LM1 || m > 0 || x < 0 || fabs(u) > 0.4){
		y = sqrt(fabs(u));
		z = sqrt(qsqfm1 + qq*xx);
		qx = q*x;
		if (qx <= 0)					{ A = z - qx;		B = q*z - x; }
		if (qx < 0 && LM1)				{ AA = qsqfm1 / A;	BB = qsqfm1*(qq*u - xx) / B; }
		if ((qx == 0 && LM1) || qx > 0)	{ AA = z + qx;		BB = q*z + x; }
		if (qx > 0)						{ A = qsqfm1 / AA;	B = qsqfm1*(qq*u - xx) / BB; }
		if (!LM1){
			if (qx*u >= 0)
				g = x*z + q*u;
			else
				g = (xx - qq*u) / (x*z - q*u);
			f = A*y;
			if (x <= 1)
				t = m*pi + atan2(f, g);
			else{
				if (f > 0.4)
					t = log(f + g);
				else{
					fg1 = f / (g + 1);
					tmp1 = 2 * fg1;
					tmp2 = 1;
					t = tmp1;
					while (1){
						t_old = t;
						tmp1 = tmp1*fg1*fg1;
						tmp2 = tmp2 + 2;
						t = t + tmp1 / tmp2;
						if (fabs(t - t_old) < 1e-8)
							break;
					}
				}
			}
			t = 2 * (t / y + B) / u;
			if (L1 && z != 0){
				dt = (3 * x*t - 4 * (A + qx*qsqfm1) / z) / u;
				if (L2)
					ddt = (3 * t + 5 * x*dt + 4 * pow(q / z, 3)*qsqfm1) / u;
				if (L3)
					dddt = (8 * dt + 7 * x*ddt - 12 * pow(q / z, 5)*x*qsqfm1) / u;
			}
		}
		else{
			dt = B;
			ddt = BB;
			dddt = AA;
		}
	}
	else{
		double u0 = 1, u1 = 1, u2 = 1, u3 = 1, l = 0, P;
		double tq, tqsum, ttmold, tterm, tqterm;
		tmp1 = 4;
		tq = q*qsqfm1;
		if (q < 0.5)
			tqsum = 1 - q*q*q;
		else
			tqsum = (q + 1 / (1 + q))*qsqfm1;
		ttmold = tmp1 / 3;
		t = ttmold*tqsum;
		while (1){
			l = l + 1;
			P = l;
			u0 = u0*u;
			if (L1 && l>1)
				u1 = u1*u;
			if (L2 && l > 2)
				u2 = u2*u;
			if (L3 && l > 3)
				u3 = u3*u;
			tmp1 = tmp1*(P - 0.5) / P;
			tq = tq*q*q;
			tqsum = tqsum + tq;
			t_old = t;
			tterm = tmp1 / (2 * P + 3);
			tqterm = tterm*tqsum;
			t = t - u0*((1.5*P + 0.25)*tqterm / (P*P - 0.25) - ttmold*tq);
			ttmold = tterm;
			tqterm = tqterm*P;
			if (L1)
				dt = dt + tqterm*u1;
			if (L2)
				ddt = ddt + tqterm*u2*(P - 1);
			if (L3)
				dddt = dddt + tqterm*u3*(P - 1)*(P - 2);
			if (l >= N && fabs(t - t_old) < 1e-8)
				break;
		}
		if (L3)
			dddt = 8 * x*(1.5*ddt - xx*dddt);
		if (L2)
			ddt = 2 * (2 * xx*ddt - dt);
		if (L1)
			dt = -2 * x*dt;
		t = t / xx;
	}

	tval[0] = t; tval[1] = dt; tval[2] = ddt; tval[3] = dddt;
}

void XLAMB(int m, double q, double qsqfm1, double tin, double *xval){
	//xval = [N,x,xp]
	double N, x, xp;
	double t0, tmin, t, dt, ddt, dddt, ddt2;
	double tval[4];
	double thr2, w, xm, xm_old;
	double tdiff, tdiffm, tdiff0;

	thr2 = atan2(qsqfm1, 2 * q) / pi;
	if (m == 0){
		N = 1;
		TLAMB(m, q, qsqfm1, 0, 0, tval);
		t0 = tval[0];
		tdiff = tin - t0;
		if (tdiff <= 0)
			x = t0*tdiff / (-4 * tin);
		else{
			x = -tdiff / (tdiff + 4);
			w = x + 1.7*sqrt(2 * (1 - thr2));
			if (w < 0)
				x = x - pow(-w, 1 / 16)*(x + sqrt(tdiff / (tdiff + 1.5*t0)));
			w = 4 / (4 + tdiff);
			x = x*(1 + x*(0.5*w - 0.03*sqrt(w)*x));
		}
	}
	else{
		xm = 1 / (1.5*(m + 0.5)*pi);
		if (thr2 < 0.5)
			xm = xm*pow(2 * thr2, 1 / 8);
		if (thr2 >= 0.5)
			xm = xm*(2 - pow(2 - 2 * thr2, 1 / 8));
		for (int ind = 1; ind <= 12; ind++){
			if (ind == 12){
				xval[0] = -1;
				return;
			}
			TLAMB(m, q, qsqfm1, xm, 3, tval);
			tmin = tval[0]; dt = tval[1]; ddt = tval[2]; dddt = tval[3];
			if (ddt == 0)
				break;
			xm_old = xm;
			xm = xm - dt*ddt / (ddt*ddt - dt*dddt / 2);
			if (fabs(xm_old / xm - 1) <= 3e-7)
				break;
		}
		tdiffm = tin - tmin;
		if (tdiffm < 0){
			xval[0] = 0;
			return;
		}
		else if (fabs(tdiffm)<1e-8){
			xval[0] = 1;
			xval[1] = xm;
			return;
		}
		else{
			N = 3;
			if (ddt == 0)
				ddt = 6 * m*pi;
			x = sqrt(tdiffm / (ddt / 2 + tdiffm / (1 - xm) / (1 - xm)));
			w = xm + x;
			w = w * 4 / (4 + tdiffm) + (1 - w)*(1 - w);
			x = x*(1 - (1 + m + 1 * (thr2 - 0.5)) / (1 + 0.15*m)*x*(0.5*w + 0.03*x*sqrt(w))) + xm;
			ddt2 = ddt / 2;
			if (x >= 1){
				N = 1;
			}
			else{
				for (int ind = 1; ind <= 3; ind++){
					TLAMB(m, q, qsqfm1, x, 2, tval);
					t = tval[0]; dt = tval[1]; ddt = tval[2]; dddt = tval[3];
					t = tin - t;
					if (dt != 0)
						x = x + t*dt / (dt*dt + t*ddt / 2);
				}
				N = 2;
				xp = x;
			}
			TLAMB(m, q, qsqfm1, 0, 0, tval);
			t0 = tval[0];
			tdiff0 = t0 - tmin;
			tdiff = tin - t0;
			if (tdiff <= 0)
				x = xm - sqrt(tdiffm / (ddt2 - tdiffm*(ddt2 / tdiff0 - 1 / xm / xm)));
			else{
				x = -tdiff / (tdiff + 4);
				w = x + 1.7*sqrt(2 * (1 - thr2));
				if (w < 0)
					x = x - pow(-w, 1 / 16)*(x + sqrt(tdiff / (tdiff + 1.5*t0)));
				w = 4 / (4 + tdiff);
				x = x*(1 + (1 + m + 0.24*(thr2 - 0.5)) / (1 + 0.15*m)*x*(0.5*w - 0.03*x*sqrt(w)));
				if (x <= -1){
					N = N - 1;
					if (N == 1)
						x = xp;
				}
			}
		}
	}
	for (int ind = 1; ind <= 3; ind++){
		TLAMB(m, q, qsqfm1, x, 2, tval);
		t = tval[0]; dt = tval[1]; ddt = tval[2];
		t = tin - t;
		if (dt != 0)
			x = x + t*dt / (dt*dt + t*ddt / 2);
	}

	xval[0] = N; xval[1] = x;
	if (N == 2)
		xval[2] = xp;
}

void VLAMB(double rm1, double rm2, double thetaf_deg, double tf, double *sol_gd){
	//sol_gd = [nsol,x,vr1,vt1,vr2,vt2,...]
	double thetaf, r1r2, cc, c, s, qsqfm1, q, rho, sig, T, GMS;
	double N, x1, x2;
	double qzminx, qzplx, zplqx, tmp;
	double tval[4], xval[3];
	int m, ind_gd;

	ind_gd = 0;
	thetaf = thetaf_deg*d2r;
	r1r2 = rm1*rm2;
	cc = rm1*rm1 + rm2*rm2 - 2 * r1r2*cos(thetaf);
	c = sqrt(cc);
	s = (rm1 + rm2 + c) / 2;
	qsqfm1 = c / s;
	q = sqrt(r1r2)*cos(thetaf / 2) / s;
	if (c != 0){
		rho = (rm1 - rm2) / c;
		sig = 4 * r1r2*sin(thetaf / 2)*sin(thetaf / 2) / cc;
	}
	else{
		rho = 0;
		sig = 1;
	}
	GMS = sqrt(mu*s / 2);
	T = 4 * GMS*tf / s / s;
	m = 0;
	while (1){
		XLAMB(m, q, qsqfm1, T, xval);
		N = xval[0];
		if (N == 0)
			break;
		else{
			x1 = xval[1];
			sol_gd[ind_gd * 5 + 1] = x1;
			TLAMB(m, q, qsqfm1, x1, -1, tval);
			qzminx = tval[1]; qzplx = tval[2]; zplqx = tval[3];
			tmp = GMS*zplqx*sqrt(sig);
			sol_gd[ind_gd * 5 + 2] = GMS*(qzminx - qzplx*rho) / rm1;
			sol_gd[ind_gd * 5 + 3] = tmp / rm1;
			sol_gd[ind_gd * 5 + 4] = -GMS*(qzminx + qzplx*rho) / rm2;
			sol_gd[ind_gd * 5 + 5] = tmp / rm2;
			ind_gd = ind_gd + 1;
			if (N == 2){
				x2 = xval[2];
				sol_gd[ind_gd * 5 + 1] = x2;
				TLAMB(m, q, qsqfm1, x2, -1, tval);
				qzminx = tval[1]; qzplx = tval[2]; zplqx = tval[3];
				tmp = GMS*zplqx*sqrt(sig);
				sol_gd[ind_gd * 5 + 2] = GMS*(qzminx - qzplx*rho) / rm1;
				sol_gd[ind_gd * 5 + 3] = tmp / rm1;
				sol_gd[ind_gd * 5 + 4] = -GMS*(qzminx + qzplx*rho) / rm2;
				sol_gd[ind_gd * 5 + 5] = tmp / rm2;
				ind_gd = ind_gd + 1;
			}
		}
		m = m + 1;
	}
	sol_gd[0] = ind_gd;
}

double NORM(double *x) {
	double res;
	res = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
	return res;
}

void NORMA(double *x, double *n) {
	double nrm;
	nrm = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
	n[0] = x[0] / nrm;
	n[1] = x[1] / nrm;
	n[2] = x[2] / nrm;
}

double DOT(double* x, double* y) {
	return(x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}

void CROSS(double *x, double *y, double *z) {
	z[0] = x[1] * y[2] - x[2] * y[1];
	z[1] = x[2] * y[0] - x[0] * y[2];
	z[2] = x[0] * y[1] - x[1] * y[0];
}

double MIN(double *x, int x_len){
	double minval = x[0];
	if (x_len>1){
		for (int i = 1; i < x_len; i++){
			if (x[i] < minval)
				minval = x[i];
		}
	}
	return minval;
}

double MAX(double *x, int x_len){
	double maxval = x[0];
	if (x_len>1){
		for (int i = 1; i < x_len; i++){
			if (x[i] > maxval)
				maxval = x[i];
		}
	}
	return maxval;
}

int issame(double *x1, double *x2, double *x3, double *x4){
	double td[4], ttr[4], J[4];
	
	td[0] = x1[0];	ttr[0] = x1[1];	J[0] = x1[2];
	td[1] = x2[0];	ttr[1] = x2[1];	J[1] = x2[2];
	td[2] = x3[0];	ttr[2] = x3[1];	J[2] = x3[2];
	td[3] = x4[0];	ttr[3] = x4[1];	J[3] = x4[2];

	if (MAX(td, 4) - MIN(td, 4) > eps_same_t)
		return 0;
	if (MAX(ttr, 4) - MIN(ttr, 4) > eps_same_t)
		return 0;
	if (MAX(J, 4) - MIN(J, 4) > eps_same_J)
		return 0;
	else
		return 1;
}

void xupdate(double *xopt, double *xshort, double *xlong){
	int ind;
	if (xopt[3] == 0){
		if (xopt[4] == 1 && xopt[2] < xshort[2]){
			for (ind = 0; ind < 3; ind++)
				xshort[ind] = xopt[ind];
		}
		else if (xopt[4] == -1 && xopt[2] < xlong[2]){
			for (ind = 0; ind < 3; ind++)
				xlong[ind] = xopt[ind];
		}
	}
}
