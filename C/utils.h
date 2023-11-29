#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void print_array( double *A, int n, int fill ) {
	
	int i;

	if ( fill == 1 ) {
		for ( i = 0; i < n; i++ ) {
				A[i] = i;
				printf("%lf  \n", A[i]);
		}
	} else {
		for ( i = 0; i < n; i++ ) {
				printf("%lf  \n", A[i]);
		}
	}
	printf("\n");
}


void print_2d_array( double **A, int num_rows, int num_cols, int fill ) {
	
	int i, j;

	if ( fill == 1 ) {
		for ( i = 0; i < num_rows; i++ ) {
			for ( j = 0; j < num_cols; j++) {
				A[i][j] = i + j;
				printf("%lf  " ,A[i][j]);
			}
			printf("\n");
		}
	} else {
		for ( i = 0; i < num_rows; i++ ) {
			for ( j = 0; j < num_cols; j++) {
				printf("%lf  " ,A[i][j]);
			}
			printf("\n");
		}
	}
	printf("\n");
}


void write_matrix( double **A, int num_rows, int num_cols, char *filename ) {

	int i, j;

	FILE *ft = fopen( filename, "w");

	for ( i = 0; i < num_rows; i++ ) {
		for (j = 0; j < num_cols; j++) {
			fprintf(ft, "%f ",  A[i][j] );
		}
		fprintf(ft, "\n");
	}

	fclose(ft);
}



void write_array( double *A, int n, char *filename ) {

	int i;

	FILE *ft = fopen( filename, "w");

	for ( i = 0; i < n; i++ ) {
		fprintf(ft, "%f ",  A[i]);
	}

	fclose(ft);
}


double dot_prod(double *x, double *y, int n) {

	double dotprod;
	int i;

	dotprod = 0.0;
	for (i = 0; i < n; i++) {
		dotprod += x[i]*y[i];
	}

	return dotprod;
}