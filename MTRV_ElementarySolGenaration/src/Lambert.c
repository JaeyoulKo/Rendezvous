#include <math.h>
#include <stdlib.h>
#include "Lambert.h"
#include "parameter.h"

double lb_ttr, ub_ttr;
double Q1[6], Q2[6], tmax;
double Jval[3], sol_gd[10000], x_local[30000];

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
