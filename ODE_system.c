
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

//compile with:
// gcc -lgsl -lm -lblas -o O ODE_system.c -L /usr/local/lib -I /usr/local/include/

int func_ode(double t, const double y[], double f[], void *params);
void simulate_ODE(double *parms, double *y0, unsigned int Tf);
void simulate_ODE_suppressed(double *parms, double *y0, unsigned int Tf);
double fitness(double *params, double R, int species);
void linspace(double min, double max, int points, double *c);

int main() {
	FILE *file = fopen("CASE0_PDE.dat","w");
	double *parms = (double*)calloc(12,sizeof(double));
	parms[0] = 1.0; //r
	parms[1] = 6.0; //K
	parms[2] = 0.5; //sigma1
	parms[3] = 0.5; //sigma2
	parms[4] = 3.0; //H1
	parms[5] = 2.5; //H2
	parms[6] = 1.0; //Imax1
	parms[7] = 1.1; //Imax2
	parms[8] = 0.015; //d1 
	parms[9] = 0.04; //d2
	parms[10] = 0.1; //T
	parms[11] = 1e-4; //sm
	
	double y0[3] = {1,1,1};
	unsigned int Tf = 10000;
	simulate_ODE_suppressed(parms, y0, Tf);
	//y0[1] = 10;
	Tf = 40000;
	simulate_ODE_suppressed(parms, y0, Tf);
	Tf = 1000;
	simulate_ODE(parms, y0, Tf);
	
	
	double Rmin1 = 0.035;//sp2 alone.
	double Rmax1 = 3.98;
	double Rmin2 = 0.09;//sp1 alone.
	double Rmax2 = 3.27;
	
	int N = 10000;
	double R_vec[N];
	linspace(0.0, 4.0, N, R_vec);
	for (int i = 0; i < N; i++) {
		if (R_vec[i] >= Rmin1 && R_vec[i] <= Rmax1) {
			if (R_vec[i] >= Rmin2 && R_vec[i] <= Rmax2) {
				fprintf(file,"%f %f %f\n",R_vec[i], fitness(parms, R_vec[i], 1), fitness(parms, R_vec[i], 2));
			} else {
				fprintf(file,"%f %f %f\n",R_vec[i], fitness(parms, R_vec[i], 1), nan(""));
			}
		} else {
			if (R_vec[i] >= Rmin2 && R_vec[i] <= Rmax2) {
				fprintf(file,"%f %f %f\n",R_vec[i], nan(""), fitness(parms, R_vec[i], 2));
			} else {
				fprintf(file,"%f %f %f\n",R_vec[i], nan(""), nan(""));
			}
		}
	}
	fclose(file);
	
	return 0;
}

void linspace(double min, double max, int points, double *c) {
	double h = (max-min)/(points-1);
	for (int i = 0; i < points; i++) {
		*(c+i) = min + h*i;
	}
}

double fitness(double *params, double R, int species) {
	double r = *((double*)params);
	double K = *((double*)params+1);
	double sigma1 = *((double*)params+2);
	double sigma2 = *((double*)params+3);
	double H1 = *((double*)params+4);
	double H2 = *((double*)params+5);
	double Imax1 = *((double*)params+6);
	double Imax2 = *((double*)params+7);
	double d1 = *((double*)params+8);
	double d2 = *((double*)params+9);
	double T = *((double*)params+10);
	double sm = *((double*)params+11);
	
	double fit = 0.0;
	if (species == 1) {
		fit = 0.05*(sigma1*R*Imax1/(H1+R)-T-d1);
	} else {
		fit = 0.05*(sigma2*R*Imax2/(H2+R)-T-d2);
	}
    return fit;
	
}


int func_ode(double t, const double y[], double f[], void *params) {
	double r = *((double*)params);
	double K = *((double*)params+1);
	double sigma1 = *((double*)params+2);
	double sigma2 = *((double*)params+3);
	double H1 = *((double*)params+4);
	double H2 = *((double*)params+5);
	double Imax1 = *((double*)params+6);
	double Imax2 = *((double*)params+7);
	double d1 = *((double*)params+8);
	double d2 = *((double*)params+9);
	double T = *((double*)params+10);
	double sm = *((double*)params+11);
	
	double R = y[0];
	double N1 = y[1];
	double N2 = y[2];
	f[0] = R*(r*(1.0-R/K)-sm*N1*Imax1/(H1+R)-sm*N2*Imax2/(H2+R));
	f[1] = N1*(sigma1*R*Imax1/(H1+R)-T-d1);
	f[2] = N2*(sigma2*R*Imax2/(H2+R)-T-d2);
    return GSL_SUCCESS;
}


void simulate_ODE(double *parms, double *y0, unsigned int Tf) {
	
    gsl_odeiv2_system sys = {func_ode, (void *)0, 3, (void *)parms};
    gsl_odeiv2_driver *d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,1e-8, 1e-8, 0.0);
	
	double t = 0.0;
    double y[3];
	double f[3];
	y[0] = y0[0];
	y[1] = y0[1];
	y[2] = y0[2];
	
	func_ode(t,y,f,(void *)parms);
    printf("%f %f %f %f\n", 0.0, y[0], y[1]*parms[11], y[2]*parms[11]);

    for (int i = 1; i <= Tf; i++)
    {
        double ti = i;
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
        
        if (status != GSL_SUCCESS)
        {
            printf ("error, return value=%d\n", status);
            break;
        }
		printf("%f %f %f %f\n", ti, y[0], y[1]*parms[11], y[2]*parms[11]);
    }
    gsl_odeiv2_driver_free (d);
}

void simulate_ODE_suppressed(double *parms, double *y0, unsigned int Tf) {
	
    gsl_odeiv2_system sys = {func_ode, (void *)0, 3, (void *)parms};
    gsl_odeiv2_driver *d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,1e-8, 1e-8, 0.0);
	
	double t = 0.0;
    double y[3];
	double f[3];
	y[0] = y0[0];
	y[1] = y0[1];
	y[2] = y0[2];
	
	func_ode(t,y,f,(void *)parms);

    for (int i = 1; i <= Tf; i++)
    {
        double ti = i;
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
        
        if (status != GSL_SUCCESS)
        {
            printf ("error, return value=%d\n", status);
            break;
        }
    }
	y0[0] = y[0];
	y0[1] = y[1];
	y0[2] = y[2];
    gsl_odeiv2_driver_free (d);
}

