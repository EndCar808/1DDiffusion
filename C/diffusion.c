#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "utils.h"	


#define PI 3.141592654
#define T_STEPS 1001
#define S_STEPS 1000


void thomas_algo(double *a, double *b, double *c, double *d, double *U, int indx);
double exact_soln(double x, double t, double t0, double alpha, char init_cond);


int main ( int argc, char **argv ) {

	/*--- VARIABLE DECLARATIONS ---*/
	double dx, dt, alpha, T, x, r, time, t0;
	double dot1, dot2;
	int t, i, j, k;
	int X_MAX;

	// X_MAX = S_STEPS + 1;
	X_MAX = S_STEPS;

	
	/* Allocate Arrays */
	// matrix arrays
	double *b1      = malloc( X_MAX * sizeof(double) );
	double *b2      = malloc( X_MAX * sizeof(double) );
	double *a       = malloc( (X_MAX) * sizeof(double) );
	double *c       = malloc( (X_MAX) * sizeof(double) );
	// first linear problem arrays
	double *y       = malloc( X_MAX * sizeof(double) );
	double *b       = malloc( X_MAX * sizeof(double) );
	// second linear problem arrays
	double *z       = malloc( X_MAX * sizeof(double) );
	double *u       = malloc( X_MAX * sizeof(double) );
	// // solution arrays
	double *U       = malloc( X_MAX * sizeof(double) );
	double *v       = malloc( X_MAX * sizeof(double) );
	double *tmp_dbl_ptr;


	double **Exact = malloc( ( T_STEPS ) * sizeof(*Exact) + ( T_STEPS ) * ( X_MAX )* sizeof(**Exact));
	tmp_dbl_ptr = Exact + ( T_STEPS );
	for (i = 0; i < ( T_STEPS ); i++, tmp_dbl_ptr += X_MAX ) {
		(Exact)[i] = tmp_dbl_ptr;
	}

	double **Soln  = malloc( ( T_STEPS ) * sizeof(*Soln) + ( T_STEPS ) * ( X_MAX )* sizeof(**Soln));
	tmp_dbl_ptr = Soln + ( T_STEPS );
	for (i = 0; i < ( T_STEPS ); i++, tmp_dbl_ptr += X_MAX ) {
		(Soln)[i] = tmp_dbl_ptr;
	}
	

	/*--- INITIALIZE VARIABLES ---*/
	T     = 1.0;
	t0    = 0.0;
	dx    = (2.0*PI)/(double)X_MAX;
	dt    = (T - t0)/(double)T_STEPS;
	time  = t0;
	alpha = 0.5;
	r     = (alpha*dt)/(2.0*dx*dx);

	
	/*--- INITIAL CONDITIONS & SETUP---*/
	for ( i = 0; i < X_MAX; i++) {
		// ICs
		x           = 0.0 + i*dx;

		// printf("%lf \n", x);
		U[i]        = exact_soln(x, t0, t0, alpha, 's');

		//printf("%lf \n", U[i]);

		Soln[0][i]  = U[i];
		Exact[0][i] = exact_soln(x, t0, t0, alpha, 's');
		
		// Fill arrays
		a[i]  = -r;
		b1[i] = 1.0 + 2.0*r;
		b2[i] = 1.0 - 2.0*r;
		c[i]  = -r;
		u[i]  = 0.0;
		v[i]  = 0.0;
		if(i == 0) {
			u[i] = -b1[i];
			v[i] = 1.0;
			b1[i] -= u[i];
		}
		if(i == (X_MAX - 1)) {
			u[i] = c[i];
			v[i] = -a[0]/b1[0];
			b1[i] -= u[i]*v[i]; 
		}
	}	
	


	/*--- INTEGRATE OVER T ---*/
	for ( t = 1; t < T_STEPS; t++) {

		// b vector
		b[0] = r*U[X_MAX - 1] + b2[0]*U[0] + r*U[1];
		for( i = 1; i < (X_MAX - 1); i++) {
			b[i] = r*U[i - 1] + b2[i]*U[i] + r*U[i + 1];
		}
		b[X_MAX - 1] = r*U[X_MAX - 2] + b2[X_MAX - 1]*U[X_MAX - 1] + r*U[0];


		thomas_algo(a, b1, c, b, y, X_MAX);
		thomas_algo(a, b1, c, u, z, X_MAX);


		//printf("u = \n");
		time = t0 + t * dt;
		dot1 = dot_prod(v, y, X_MAX);
		dot2 = 1.0 + dot_prod(v, z, X_MAX);
		for (i = 0; i < X_MAX; i++ ) {
			
			U[i] = y[i] - (dot1/dot2)*z[i];
			//printf("%lf\n", U[i]);
			// save answers
			x           = 0.0 + i*dx;
			Soln[t][i]  = U[i];
			Exact[t][i] = exact_soln(x, time, t0, alpha, 's');

		}
		//printf("\n\n");
	}

	FILE *v1 = fopen( "../Data/C-Data/variables.txt", "w");
	fprintf(v1, "%d\n", X_MAX);
	fprintf(v1, "%lf\n", dx);
	fprintf(v1, "%lf\n", dt);
	fprintf(v1, "%lf\n", alpha);
	fprintf(v1, "%lf\n", T);
	fclose(v1);

	write_matrix(Soln, T_STEPS, X_MAX, "../Data/C-Data/solution.txt");
	write_matrix(Exact, T_STEPS, X_MAX, "../Data/C-Data/exact.txt");

	//print_array( b, X_MAX, 0);
	//print_2d_array(A, X_MAX, X_MAX, 0);


	free(b1);


	return 0;
}


void thomas_algo(double *a, double *b, double *c, double *d, double *U, int indx) {

	int i;
	double tmp;

	double *cc = malloc((indx - 1)*sizeof(*cc));
	double *dd = malloc(indx*sizeof(*dd));

	// forward sweep
	cc[0] = c[0]/b[0];
	dd[0] = d[0]/b[0];
	for (i = 1; i < indx; ++i) {
		tmp = 1.0 / (b[i] - a[i]*cc[i - 1]); // recall a has length index - 1 so we shift back by 1
		if (i < (indx - 1)) {
			cc[i] = c[i]*tmp;
		}
		dd[i] = (d[i] - a[i]*dd[i - 1])*tmp;
	}

	// back substitution
	U[indx - 1] = dd[indx - 1];
	for (i = indx - 2; i >= 0; i--) {
		U[i] = dd[i] - cc[i]*U[i + 1];
	}	

	free(cc);
	free(dd);
}


double exact_soln(double x, double t, double t0, double alpha, char init_cond) {

	double soln;


	switch(init_cond) {
		case 'p' :
			soln = sqrt(t0 / t)*exp(-(x*x)/(4.0*alpha*t));
			break;
		case 's' :
			soln = sin(x)*exp(-alpha*t);
			break;
		case 'c':
			soln = cos(x)*exp(-alpha*t);
			break;
		default:
			soln = sin(x)*exp(-alpha*t);
	}

	return soln;
}