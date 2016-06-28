#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "math_lib.h"

void matsavefile(char *str, double **M, int m, int n)
{
	int i,j;
	FILE *f;
	f=fopen( str, "w+");

	for( i=0; i<m; i++ ){
		for( j=0; j<n; j++ ){
			fprintf(f, "%09.7f ", M[i][j] );
		}
		fprintf(f, "\n" );
	}
	fclose(f);
}


void matsavefilex(char *str, double **M, int ms, int mf, int ns, int nf)
{
	int i,j;
	FILE *f;
	f=fopen(str, "w+");

	for( i=ms; i<=mf; i++ ){
		for( j=ns; j<=nf; j++ ){
			fprintf(f, "%09.7f ", M[i][j] );
		}
		fprintf(f, "\n" );
	}
	fclose(f);
}


double **smatrix(unsigned int m, unsigned int n)
{
	unsigned int i;
	double **M;

	M = (double **) malloc( m*sizeof(double *) );
	M[0] = (double *) malloc((m*n)*sizeof(double));

	for(i=1; i<m; i++) M[i]=M[i-1]+n;
	return M;
}

double *svector(unsigned int m)
{
	double *V;

	V=(double *)malloc( m*sizeof(double));
	return V;
}


unsigned int matsvdinvall( double **mat, int m, int n, double singular_tollerance, double *det, double **out_mat)
{
	double **mat_in;
	double *S;
	double **V;
	int i, j, k;
	int dof, rank;
	unsigned int *index_dof;

	dof = (m>n?n:m);
	S=dvector(1,n);
	V=dmatrix(1,n,1,n);
	index_dof = (unsigned int *)malloc(dof*sizeof(unsigned int));

	// copy input matrix
	mat_in = dmatrix(1, m, 1, n);
	for(i=0; i<m; i++ )
		for( j=0; j<n; j++)
			mat_in[i+1][j+1] = mat[i][j];

	// SVD
	svdcmp(mat_in, m, n, S, V );
	rank = 0;
	*det = 1;
	for( i=0; i<n ; i++ ){
		if( fabs(S[i+1]) > singular_tollerance ){
			index_dof[rank++] = i;
			*det *= S[i+1];
		}

	}

	// inverse SVD
	for( j=0; j<n; j++ ){
		for( k=0; k<rank; k++ ){
			V[j+1][ index_dof[k]+1 ] /=S[ index_dof[k] +1];
		}
	}

	for( i=0; i<n; i++ ){
		for( j=0; j<m; j++ ){
			out_mat[i][j] = 0;
			for( k=0; k<rank; k++ ){
				out_mat[i][j] += V[i+1][ index_dof[k]+1 ]*mat_in[j+1][ index_dof[k]+1 ];
			}
		}
	}

        free(index_dof);
	return rank;
}

unsigned int matsvdinv( double **mat, int m, int n, double **out_mat)
{
	double **mat_in;
	double *S;
	double **V;
	int i, j, k;
	int dof, rank;
	unsigned int *index_dof;

	dof = (m>n?n:m);
	S=dvector(1,n);
	V=dmatrix(1,n,1,n);
	index_dof = (unsigned int *)malloc(dof*sizeof(unsigned int));

	// copy input matrix
	mat_in = dmatrix(1, m, 1, n);
	for(i=0; i<m; i++ )
		for( j=0; j<n; j++)
			mat_in[i+1][j+1] = mat[i][j];

	// SVD
	svdcmp(mat_in, m, n, S, V );
	rank = 0;
	for( i=0; i<n ; i++ ){
		if( fabs(S[i+1]) > SINGULAR_TOLERANCE ) index_dof[rank++] = i;
	}

	// inverse SVD
	for( j=0; j<n; j++ ){
		for( k=0; k<rank; k++ ){
			V[j+1][ index_dof[k]+1 ] /=S[ index_dof[k] +1];
		}
	}

	for( i=0; i<n; i++ ){
		for( j=0; j<m; j++ ){
			out_mat[i][j] = 0;
			for( k=0; k<rank; k++ ){
				out_mat[i][j] += V[i+1][ index_dof[k]+1 ]*mat_in[j+1][ index_dof[k]+1 ];
			}
		}
	}

        free(index_dof);
	return rank;
}

unsigned int matsvdinvDLS( double **mat, int m, int n, double lamda, double **out_mat)
{
	double **mat_in;
	double *S;
	double **V;
	int i, j, k;
	int dof, rank;
	unsigned int *index_dof;

	dof = (m>n?n:m);
	S=dvector(1,n);
	V=dmatrix(1,n,1,n);
	index_dof = (unsigned int *)malloc(dof*sizeof(unsigned int));

	// copy input matrix to mat_in
	mat_in = dmatrix(1, m, 1, n);
	for(i=0; i<m; i++ )
		for( j=0; j<n; j++)
			mat_in[i+1][j+1] = mat[i][j];

	// SVD
	svdcmp(mat_in, m, n, S, V );
	rank = 0;
	for( i=0; i<n ; i++ ){
		if( fabs(S[i+1]) > SINGULAR_TOLERANCE ) index_dof[rank++] = i;
	}

	// inverse SVD Dampled Least Mean Square
	for( j=0; j<n; j++ ){
		for( k=0; k<rank; k++ ){
			V[j+1][ index_dof[k]+1 ] *=( S[ index_dof[k] +1]/ (S[ index_dof[k] +1]*S[ index_dof[k] +1] +lamda*lamda) );
		}
	}

	for( i=0; i<n; i++ ){
		for( j=0; j<m; j++ ){
			out_mat[i][j] = 0;
			for( k=0; k<rank; k++ ){
				out_mat[i][j] += V[i+1][ index_dof[k]+1 ]*mat_in[j+1][ index_dof[k]+1 ];
			}
		}
	}

        free(index_dof);
	return rank;
}

void mattranspose( double **mat, int m, int n, double **out_mat)
{
    int i,j;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[j][i] = mat[i][j];
}

void matmul_vec(double **mat, int m, int n, double *rvec, double *out_vec )
{
    int i,j;

	for(i=0; i<m; i++){
		out_vec[i] = 0.0;
		for(j=0; j<n; j++)
			out_vec[i] += mat[i][j]*rvec[j];
	}
}



void matzero_c7( int m, double out_mat[][7] )
{
    int n = 7;
    int i,j;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] = 0.;
}

void matzero_c4( int m, double out_mat[][4] )
{
    int n = 4;
    int i,j;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] = 0.;
}


void matzero_c6( int m, double out_mat[][6] )
{
    int n = 6;
    int i,j;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] = 0.;
}

void matzero_c1( int m, double *out_vec )
{
    int i;

    for(i=0; i<m; i++)
		out_vec[i] = 0.;
}

void matmul_cc( int lm, int ln, int rn,  double **lmat, double **rmat, double **out_mat)
{
	int i,j,k;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
			out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c42( int lm, double lmat[][4], double rmat[][2], double out_mat[][2])
{
    int i,j,k;
    int rn = 2;
    int ln = 4;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c44( int lm, double lmat[][4], double rmat[][4], double out_mat[][4] )
{
    int i,j,k;
    int rn = 4;
    int ln = 4;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c33( int lm, double lmat[][3], double rmat[][3], double out_mat[][3] )
{
    int i,j,k;
    int rn = 3;
    int ln = 3;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c34( int lm, double lmat[][3], double rmat[][4], double out_mat[][4])
{
    int i,j,k;
    int rn = 4;
    int ln = 3;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c36( int lm, double lmat[][3], double rmat[][6], double out_mat[][6])
{
    int i,j,k;
    int rn = 6;
    int ln = 3;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c63( int lm, double lmat[][6], double rmat[][3], double out_mat[][3])
{
    int i,j,k;
    int rn = 3;
    int ln = 6;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c613( int lm, double lmat[][6], double rmat[][13], double out_mat[][13])
{
    int i,j,k;
    int rn = 13;
    int ln = 6;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c164( int lm, double lmat[][16], double rmat[][4], double out_mat[][4])
{
    int i,j,k;
    int rn = 4;
    int ln = 16;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c153( int lm, double lmat[][15], double rmat[][3], double out_mat[][3])
{
    int i,j,k;
    int rn = 3;
    int ln = 15;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c133( int lm, double lmat[][13], double rmat[][3], double out_mat[][3])
{
    int i,j,k;
    int rn = 3;
    int ln = 13;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c316( int lm, double lmat[][3], double rmat[][16], double out_mat[][16])
{
    int i,j,k;
    int rn = 16;
    int ln = 3;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c313( int lm, double lmat[][3], double rmat[][13], double out_mat[][13])
{
    int i,j,k;
    int rn = 13;
    int ln = 3;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c1515( int lm, double lmat[][15], double rmat[][15], double out_mat[][15])
{
    int i,j,k;
    int rn = 15;
    int ln = 15;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c1616( int lm, double lmat[][16], double rmat[][16], double out_mat[][16])
{
    int i,j,k;
    int rn = 16;
    int ln = 16;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c1313( int lm, double lmat[][13], double rmat[][13], double out_mat[][13])
{
    int i,j,k;
    int rn = 13;
    int ln = 13;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c66( int lm, double lmat[][6], double rmat[][6], double out_mat[][6])
{
    int i,j,k;
    int rn = 6;
    int ln = 6;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}

}

void matmul_c37( int lm, double lmat[][3], double rmat[][7], double out_mat[][7])
{
    int i,j,k;
    int rn = 7;
    int ln = 3;

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matmul_c21( int lm, double lmat[][2], double *rvec, double *out_vec)
{
    int i,k;
    int ln = 2;

    for(i=0; i<lm; i++)
    {
	out_vec[i] = 0.;
	for(k=0; k<ln; k++)
	    out_vec[i] += lmat[i][k]*rvec[k];
    }
}

void matmul_c31( int lm, double lmat[][3], double *rvec, double *out_vec)
{
    int i,k;
    int ln = 3;

    for(i=0; i<lm; i++)
    {
	out_vec[i] = 0.;
	for(k=0; k<ln; k++)
	    out_vec[i] += lmat[i][k]*rvec[k];
    }
}

void matmul_c41( int lm, double lmat[][4], double *rvec, double *out_vec)
{
    int i,k;
    int ln = 4;

    for(i=0; i<lm; i++)
    {
	out_vec[i] = 0.;
	for(k=0; k<ln; k++)
	    out_vec[i] += lmat[i][k]*rvec[k];
    }
}

void matmul_c( int row, int col,  double **lmat, double *rvec, double *out_vec)
{
    int i,k;

    for(i=0; i<row; i++)
    {
		out_vec[i] = 0.;
		for(k=0; k<col; k++)
			out_vec[i] += lmat[i][k]*rvec[k];
    }
}



void matmul_c61( int lm, double lmat[][6], double *rvec, double *out_vec)
{
    int i,k;
    int ln = 6;

    for(i=0; i<lm; i++)
    {
	out_vec[i] = 0.;
	for(k=0; k<ln; k++)
	    out_vec[i] += lmat[i][k]*rvec[k];
    }
}

void matmul_c101( int lm, double lmat[][10], double *rvec, double *out_vec)
{
    int i,k;
    int ln = 10;

    for(i=0; i<lm; i++)
    {
	out_vec[i] = 0.;
	for(k=0; k<ln; k++)
	    out_vec[i] += lmat[i][k]*rvec[k];
    }
}

void matmul_c131( int lm, double lmat[][13], double *rvec, double *out_vec)
{
    int i,k;
    int ln = 13;

    for(i=0; i<lm; i++)
    {
	out_vec[i] = 0.;
	for(k=0; k<ln; k++)
	    out_vec[i] += lmat[i][k]*rvec[k];
    }
}

void matmul_c161( int lm, double lmat[][16], double *rvec, double *out_vec)
{
    int i,k;
    int ln = 16;

    for(i=0; i<lm; i++)
    {
	out_vec[i] = 0.;
	for(k=0; k<ln; k++)
	    out_vec[i] += lmat[i][k]*rvec[k];
    }
}

void matmul_c71( int lm, double lmat[][7], double *rvec, double *out_vec)
{
    int i,k;
    int ln = 7;

    for(i=0; i<lm; i++)
    {
	out_vec[i] = 0.;
	for(k=0; k<ln; k++)
	    out_vec[i] += lmat[i][k]*rvec[k];
    }
}

void matrmul_r33c33( double rmat[][3], double out_mat[][3])
{
    double lmat[3][3];
    int rn = 3;
    int lm = 3;
    int ln = 3;
    int i,j,k;

    matcp_c3(3, out_mat, lmat);

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matrmul_r44c44( double rmat[][4], double out_mat[][4] )
{
    double lmat[4][4];
    int rn = 4;
    int lm = 4;
    int ln = 4;
    int i,j,k;

    matcp_c4(4, out_mat, lmat);

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matrmul_r66c66( double rmat[][6], double out_mat[][6])
{
    double lmat[6][6];
    int rn = 6;
    int lm = 6;
    int ln = 6;
    int i,j,k;

    matcp_c6(6, out_mat, lmat);

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matlmul_r33c36( double lmat[][3], double out_mat[][6])
{
    double rmat[3][6];
    int lm = 3;
    int rn = 6;
    int ln = 3;
    int i, j, k;

    matcp_c6(3, out_mat, rmat);

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matlmul_r44c44( double lmat[][4], double out_mat[][4] )
{
    double rmat[4][4];
    int lm = 4;
    int rn = 4;
    int ln = 4;
    int i, j, k;

    matcp_c4(4, out_mat, rmat);

    for(i=0; i<lm; i++)
	for(j=0; j<rn; j++)
	{
	    out_mat[i][j] = 0.;
	    for(k=0; k<ln; k++)
		out_mat[i][j] += lmat[i][k]*rmat[k][j];
	}
}

void matadd_c4( int m, double mat[][4], double out_mat[][4] )
{
    int n = 4;
    int i, j;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] += mat[i][j];
}

void matadd_c3( int m, double mat[][3], double out_mat[][3])
{
    int n = 3;
    int i, j;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] += mat[i][j];
}

void matadd_c7( int m, double mat[][7], double out_mat[][7] )
{
    int n = 7;
    int i, j;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] += mat[i][j];
}

void matadd_c6( int m, double mat[][6], double out_mat[][6] )
{
    int n = 6;
    int i, j;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] += mat[i][j];
}

void matadd_c1( int m, double *vec, double *out_vec )
{
    int i;

    for(i=0; i<m; i++)
	out_vec[i] += vec[i];
}

void matadd_c( int m, double *lvec, double *rvec, double *out_vec )
{
    int i;

    for(i=0; i<m; i++)
	out_vec[i] = lvec[i] + rvec[i];
}


void matsub_c1( int m, const double *lvec, const double *rvec, double *out_vec )
{
    int i;

    for(i=0; i<m; i++)
	out_vec[i] = lvec[i] - rvec[i];
}

void matsmul_c7( int m, double scalar, double out_mat[][7] )
{
    int n=7;
    int i,j;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] *= scalar;
}

void matsmul_c4( int m, double scalar, double out_mat[][4] )
{
    int n=4;
    int i,j;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] *= scalar;
}

void matsmul_c6( int m, double scalar, double out_mat[][6] )
{
    int n=6;
    int i,j;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] *= scalar;
}

void matsmul_c1( int m, double scalar, double *out_vec )
{
    int i;

    for( i=0; i<m; i++)
	out_vec[i] *= scalar;
}

void matI_c3( double out_mat[][3] )
{
    int i,j;
    int dim = 3;

    for(i=0; i<dim; i++)
	for(j=0; j<dim; j++)
	{
	    if( i==j )
		out_mat[i][j] = 1.;
	    else
		out_mat[i][j] = 0.;
	}
}

void matI_c15( double out_mat[][15])
{
    int i,j;
    int dim = 15;

    for(i=0; i<dim; i++)
	for(j=0; j<dim; j++)
	{
	    if( i==j )
		out_mat[i][j] = 1.;
	    else
		out_mat[i][j] = 0.;
	}
}

void matI_c10( double out_mat[][10] )
{
    int i,j;
    int dim = 10;

    for(i=0; i<dim; i++)
	for(j=0; j<dim; j++)
	{
	    if( i==j )
		out_mat[i][j] = 1.;
	    else
		out_mat[i][j] = 0.;
	}
}

void matI_c16( double out_mat[][16])
{
    int i,j;
    int dim = 16;

    for(i=0; i<dim; i++)
	for(j=0; j<dim; j++)
	{
	    if( i==j )
		out_mat[i][j] = 1.;
	    else
		out_mat[i][j] = 0.;
	}
}

void matI_c13( double out_mat[][13])
{
    int i,j;
    int dim = 13;

    for(i=0; i<dim; i++)
	for(j=0; j<dim; j++)
	{
	    if( i==j )
		out_mat[i][j] = 1.;
	    else
		out_mat[i][j] = 0.;
	}
}

void matI_c6( double out_mat[][6] )
{
    int i,j;
    int dim = 6;

    for(i=0; i<dim; i++)
	for(j=0; j<dim; j++)
	{
	    if( i==j )
		out_mat[i][j] = 1.;
	    else
		out_mat[i][j] = 0.;
	}
}

void matI_c4( double out_mat[][4] )
{
    int i,j;
    int dim = 4;

    for(i=0; i<dim; i++)
	for(j=0; j<dim; j++)
	{
	    if( i==j )
		out_mat[i][j] = 1.;
	    else
		out_mat[i][j] = 0.;
	}
}

void matcp_c15( int m, double mat[][15], double out_mat[][15] )
{
    int i,j;
    int n = 15;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] = mat[i][j];
}

void matcp_c16( int m, double mat[][16], double out_mat[][16] )
{
    int i,j;
    int n = 16;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] = mat[i][j];
}

void matcp_c13( int m, double mat[][13], double out_mat[][13] )
{
    int i,j;
    int n = 13;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] = mat[i][j];
}

void matcp_c10( int m, double mat[][10], double out_mat[][10] )
{
    int i,j;
    int n = 10;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] = mat[i][j];
}

void matcp_c7( int m, double mat[][7], double out_mat[][7] )
{
    int i,j;
    int n = 7;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] = mat[i][j];
}

void matcp_c6( int m, double mat[][6], double out_mat[][6] )
{
    int i,j;
    int n = 6;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] = mat[i][j];
}

void matcp_c3( int m, double mat[][3], double out_mat[][3] )
{
    int i,j;
    int n = 3;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] = mat[i][j];
}


void matcp_c4( int m, double mat[][4], double out_mat[][4] )
{
    int i,j;
    int n = 4;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[i][j] = mat[i][j];
}

void matcp_c1( int m, double *vec, double *out_vec)
{
    int i;

    for(i=0; i<m; i++)
	out_vec[i] = vec[i];
}

void mattransp_r3c6( double mat[][6], double out_mat[][3] )
{
    int m = 3;
    int n = 6;

    int i,j;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[j][i] = mat[i][j];
}

void mattransp_r3c15( double mat[][15], double out_mat[][3] )
{
    int m = 3;
    int n = 15;

    int i,j;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[j][i] = mat[i][j];
}

void mattransp_r4c16( double mat[][16], double out_mat[][4] )
{
    int m = 4;
    int n = 16;

    int i,j;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[j][i] = mat[i][j];
}

void mattransp_r3c13( double mat[][13], double out_mat[][3] )
{
    int m = 3;
    int n = 13;

    int i,j;

    for(i=0; i<m; i++)
	for(j=0; j<n; j++)
	    out_mat[j][i] = mat[i][j];
}

void matpinv_r3c6( double mat[][6], double weight[][6], double out_mat[][3])
{
    double invW[6][6];
    double JinvW[3][6];
    double JT[6][3];
    double JinvWJT[3][3];
    double inv33[3][3];
    double invWJT[6][3];

    matinv_c6( weight, invW );
    matmul_c66(3, mat, invW, JinvW );
    mattransp_r3c6( mat, JT );
    matmul_c63( 3, JinvW, JT, JinvWJT );
    matinv_c3( JinvWJT, inv33 );
    matmul_c63( 6, invW, JT, invWJT );
    matmul_c33( 6, invWJT, inv33, out_mat );
}

void matpinv_r3c15( double mat[][15], double weight[][15], double out_mat[][3])
{
    double invW[15][15];
    double JinvW[3][15];
    double JT[15][3];
    double JinvWJT[3][3];
    double inv33[3][3];
    double invWJT[15][3];

    matinv_c15( weight, invW );
    matmul_c1515(3, mat, invW, JinvW );
    mattransp_r3c15( mat, JT );
    matmul_c153( 3, JinvW, JT, JinvWJT );
    matinv_c3( JinvWJT, inv33 );
    matmul_c153( 15, invW, JT, invWJT );
    matmul_c33( 15, invWJT, inv33, out_mat );
}

void matpinv_r3c13( double mat[][13], double weight[][13], double out_mat[][3])
{
    double invW[13][13];
    double JinvW[3][13];
    double JT[13][3];
    double JinvWJT[3][3];
    double inv33[3][3];
    double invWJT[13][3];

    matinv_c13( weight, invW );
    matmul_c1313(3, mat, invW, JinvW );
    mattransp_r3c13( mat, JT );
    matmul_c133( 3, JinvW, JT, JinvWJT );
    matinv_c3( JinvWJT, inv33 );
    matmul_c133( 13, invW, JT, invWJT );
    matmul_c33( 13, invWJT, inv33, out_mat );
}

void matpinv_r4c16( double mat[][16], double weight[][16], double out_mat[][4])
{
    double invW[16][16];
    double JinvW[4][16];
    double JT[16][4];
    double JinvWJT[4][4];
    double inv44[4][4];
    double invWJT[16][4];

    matinv_c16( weight, invW );
    matmul_c1616(4, mat, invW, JinvW );
    mattransp_r4c16( mat, JT );
    matmul_c164( 4, JinvW, JT, JinvWJT );
    matinv_c4( JinvWJT, inv44 );
    matmul_c164( 16, invW, JT, invWJT );
    matmul_c44( 16, invWJT, inv44, out_mat );
}

void matpinv_r1c13( double* row_vec, double weight[][13], double* out_col_vec)
{
    double invW[13][13];
    double JinvW[13];	    // row vector
    double JT[13];	    // column vector
    double JinvWJT;
    double inv;
    double invWJT[13];	    // column vector
    int dim = 13;
    int i, j;

    matinv_c13( weight, invW );

    for(i=0; i<dim; i++)
    {
	JinvW[i] = 0.;
	for(j=0; j<dim; j++)
	    JinvW[i] += row_vec[j]*invW[j][i];
    }

    matcp_c1( dim, row_vec, JT );

    JinvWJT = 0.;
    for(i=0; i<dim; i++)
	JinvWJT += JinvW[i] * JT[i];

    inv = 1.0/JinvWJT;

    matmul_c131( 13, invW, JT, invWJT );

    for(i=0; i<dim; i++)
	out_col_vec[i] = invWJT[i]*inv;
}

void matinv_c2( double mat[][2], double out_mat[][2] )
{
    double det;

    det = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];

    out_mat[0][0] = 1./det * mat[1][1];
    out_mat[1][0] = 1./det * (-mat[1][0]);
    out_mat[0][1] = 1./det * (-mat[0][1]);
    out_mat[1][1] = 1./det * mat[0][0];
}

void matinv_c10( double mat[][10], double out_mat[][10] )
{
    int i, j, k;
    double a1, a2;
    int m = 10;
    int n = 10;

    double temp[10][10];
    double crhs[10][10];

    matcp_c10( 10, mat, crhs );
    matI_c10( temp );
    for(k=0; k<m; k++)
    {
	int indx = matpivot_c10( k, crhs );

	// If indx = -1, Singular matrix inversion error !!

	if (indx != 0)
	    matexchrow_c10(k, indx, temp);

	a1 = crhs[k][k];
	for (j=0; j<n; j++)
	{
	    crhs[k][j] /= a1;
	    temp[k][j] /= a1;
	}
	for (i=0; i<m; i++)
	{
	    if(i!=k)
	    {
		a2 = crhs[i][k];
		for(j=0; j<n; j++)
		{
		    crhs[i][j] -= a2*crhs[k][j];
		    temp[i][j] -= a2*temp[k][j];
		}
	    }
	}
    }

    matcp_c10( 10, temp, out_mat );
}

void matinv_c15( double mat[][15], double out_mat[][15] )
{
    int i, j, k;
    double a1, a2;
    int m = 15;
    int n = 15;

    double temp[15][15];
    double crhs[15][15];

    matcp_c15( 15, mat, crhs );
    matI_c15( temp );
    for(k=0; k<m; k++)
    {
	int indx = matpivot_c15( k, crhs );

	// If indx = -1, Singular matrix inversion error !!

	if (indx != 0)
	    matexchrow_c15(k, indx, temp);

	a1 = crhs[k][k];
	for (j=0; j<n; j++)
	{
	    crhs[k][j] /= a1;
	    temp[k][j] /= a1;
	}
	for (i=0; i<m; i++)
	{
	    if(i!=k)
	    {
		a2 = crhs[i][k];
		for(j=0; j<n; j++)
		{
		    crhs[i][j] -= a2*crhs[k][j];
		    temp[i][j] -= a2*temp[k][j];
		}
	    }
	}
    }

    matcp_c15( 15, temp, out_mat );
}

void matinv_c16( double mat[][16], double out_mat[][16] )
{
    int i, j, k;
    double a1, a2;
    int m = 16;
    int n = 16;

    double temp[16][16];
    double crhs[16][16];

    matcp_c16( 16, mat, crhs );
    matI_c16( temp );
    for(k=0; k<m; k++)
    {
	int indx = matpivot_c16( k, crhs );

	// If indx = -1, Singular matrix inversion error !!

	if (indx != 0)
	    matexchrow_c16(k, indx, temp);

	a1 = crhs[k][k];
	for (j=0; j<n; j++)
	{
	    crhs[k][j] /= a1;
	    temp[k][j] /= a1;
	}
	for (i=0; i<m; i++)
	{
	    if(i!=k)
	    {
		a2 = crhs[i][k];
		for(j=0; j<n; j++)
		{
		    crhs[i][j] -= a2*crhs[k][j];
		    temp[i][j] -= a2*temp[k][j];
		}
	    }
	}
    }

    matcp_c16( 16, temp, out_mat );
}

void matinv_c13( double mat[][13], double out_mat[][13] )
{
    int i, j, k;
    double a1, a2;
    int m = 13;
    int n = 13;

    double temp[13][13];
    double crhs[13][13];

    matcp_c13( 13, mat, crhs );
    matI_c13( temp );
    for(k=0; k<m; k++)
    {
	int indx = matpivot_c13( k, crhs );

	// If indx = -1, Singular matrix inversion error !!

	if (indx != 0)
	    matexchrow_c13(k, indx, temp);

	a1 = crhs[k][k];
	for (j=0; j<n; j++)
	{
	    crhs[k][j] /= a1;
	    temp[k][j] /= a1;
	}
	for (i=0; i<m; i++)
	{
	    if(i!=k)
	    {
		a2 = crhs[i][k];
		for(j=0; j<n; j++)
		{
		    crhs[i][j] -= a2*crhs[k][j];
		    temp[i][j] -= a2*temp[k][j];
		}
	    }
	}
    }

    matcp_c13( 13, temp, out_mat );
}

void matinv_c6( double mat[][6], double out_mat[][6] )
{
    int i, j, k;
    double a1, a2;
    int m = 6;
    int n = 6;

    double temp[6][6];
    double crhs[6][6];

    matcp_c6( 6, mat, crhs );
    matI_c6( temp );
    for(k=0; k<m; k++)
    {
	int indx = matpivot_c6( k, crhs );

	// If indx = -1, Singular matrix inversion error !!

	if (indx != 0)
	    matexchrow_c6(k, indx, temp);

	a1 = crhs[k][k];
	for (j=0; j<n; j++)
	{
	    crhs[k][j] /= a1;
	    temp[k][j] /= a1;
	}
	for (i=0; i<m; i++)
	{
	    if(i!=k)
	    {
		a2 = crhs[i][k];
		for(j=0; j<n; j++)
		{
		    crhs[i][j] -= a2*crhs[k][j];
		    temp[i][j] -= a2*temp[k][j];
		}
	    }
	}
    }

    matcp_c6( 6, temp, out_mat );
}

void matinv_c4( double mat[][4], double out_mat[][4] )
{
    int i, j, k;
    double a1, a2;
    int m = 4;
    int n = 4;

    double temp[4][4];
    double crhs[4][4];

    matcp_c4( 4, mat, crhs );
    matI_c4( temp );
    for(k=0; k<m; k++)
    {
		int indx = matpivot_c4( k, crhs );

		// If indx = -1, Singular matrix inversion error !!

		if (indx != 0)
		{
			matexchrow_c4(k, indx, temp);
		}
		a1 = crhs[k][k];
		for (j=0; j<n; j++)
		{
			crhs[k][j] /= a1;
			temp[k][j] /= a1;
		}
		for (i=0; i<m; i++)
		{
			if(i!=k)
			{
			a2 = crhs[i][k];
			for(j=0; j<n; j++)
			{
				crhs[i][j] -= a2*crhs[k][j];
				temp[i][j] -= a2*temp[k][j];
			}
			}
		}
    }

    matcp_c4( 4, temp, out_mat );
}

void matinv_c3( double mat[][3], double out_mat[][3] )
{
    int i, j, k;
    double a1, a2;
    int m = 3;
    int n = 3;

    double temp[3][3];
    double crhs[3][3];

    matcp_c3( 3, mat, crhs );
    matI_c3( temp );
    for(k=0; k<m; k++)
    {
	int indx = matpivot_c3( k, crhs );

	// If indx = -1, Singular matrix inversion error !!

	if (indx != 0)
	{
	    matexchrow_c3(k, indx, temp);
	}
	a1 = crhs[k][k];
	for (j=0; j<n; j++)
	{
	    crhs[k][j] /= a1;
	    temp[k][j] /= a1;
	}
	for (i=0; i<m; i++)
	{
	    if(i!=k)
	    {
		a2 = crhs[i][k];
		for(j=0; j<n; j++)
		{
		    crhs[i][j] -= a2*crhs[k][j];
		    temp[i][j] -= a2*temp[k][j];
		}
	    }
	}
    }

    matcp_c3( 3, temp, out_mat );
}

void matnull( int m, int n, double **J, double **Ji, double **N)
{
	int i,j;
	matmul_cc(n, m, n, Ji, J, N);
	for( i=0; i<n; i++ )
	for( j=0; j<n; j++ )
		N[i][j] = -N[i][j];

	for( i=0; i<n; i++ )
		N[i][i] = 1.0+N[i][i];
}

int matpivot_c15( int nm, double mat[][15] )
{
    int k = nm;
    double amax, temp;
    int i;
    int m = 15;

    amax = -1.;
    for (i=nm; i<m; i++)
	if ( (temp = dabs(mat[i][nm])) > amax && temp != 0.0)
	{
	    amax = temp;
	    k = i;
	}
    if (mat[k][nm] == (double)(0))
	return -1;
    if (k != nm)
    {
	matexchrow_c15( k, nm, mat );
	return k;
    }
    return 0;
}

int matpivot_c16( int nm, double mat[][16] )
{
    int k = nm;
    double amax, temp;
    int i;
    int m = 16;

    amax = -1.;
    for (i=nm; i<m; i++)
	if ( (temp = dabs(mat[i][nm])) > amax && temp != 0.0)
	{
	    amax = temp;
	    k = i;
	}
    if (mat[k][nm] == (double)(0))
	return -1;
    if (k != nm)
    {
	matexchrow_c16( k, nm, mat );
	return k;
    }
    return 0;
}

int matpivot_c13( int nm, double mat[][13] )
{
    int k = nm;
    double amax, temp;
    int i;
    int m = 13;

    amax = -1.;
    for (i=nm; i<m; i++)
	if ( (temp = dabs(mat[i][nm])) > amax && temp != 0.0)
	{
	    amax = temp;
	    k = i;
	}
    if (mat[k][nm] == (double)(0))
	return -1;
    if (k != nm)
    {
	matexchrow_c13( k, nm, mat );
	return k;
    }
    return 0;
}

int matpivot_c10( int nm, double mat[][10] )
{
    int k = nm;
    double amax, temp;
    int i;
    int m = 10;

    amax = -1.;
    for (i=nm; i<m; i++)
	if ( (temp = dabs(mat[i][nm])) > amax && temp != 0.0)
	{
	    amax = temp;
	    k = i;
	}
    if (mat[k][nm] == (double)(0))
	return -1;
    if (k != nm)
    {
	matexchrow_c10( k, nm, mat );
	return k;
    }
    return 0;
}

int matpivot_c6( int nm, double mat[][6] )
{
    int k = nm;
    double amax, temp;
    int i;
    int m = 6;

    amax = -1.;
    for (i=nm; i<m; i++)
	if ( (temp = dabs(mat[i][nm])) > amax && temp != 0.0)
	{
	    amax = temp;
	    k = i;
	}
    if (mat[k][nm] == (double)(0))
	return -1;
    if (k != nm)
    {
	matexchrow_c6( k, nm, mat );
	return k;
    }
    return 0;
}

int matpivot_c4( int nm, double mat[][4] )
{
    int k = nm;
    double amax, temp;
    int i;
    int m = 4;

    amax = -1.;
    for (i=nm; i<m; i++)
	if ( (temp = dabs(mat[i][nm])) > amax && temp != 0.0)
	{
	    amax = temp;
	    k = i;
	}
    if (mat[k][nm] == (double)(0))
	return -1;
    if (k != nm)
    {
	matexchrow_c4( k, nm, mat );
	return k;
    }
    return 0;
}

int matpivot_c3( int nm, double mat[][3] )
{
    int k = nm;
    double amax, temp;
    int i;
    int m = 3;

    amax = -1.;
    for (i=nm; i<m; i++)
	if ( (temp = dabs(mat[i][nm])) > amax && temp != 0.0)
	{
	    amax = temp;
	    k = i;
	}
    if (mat[k][nm] == (double)(0))
	return -1;
    if (k != nm)
    {
	matexchrow_c3( k, nm, mat );
	return k;
    }
    return 0;
}

void matexchrow_c15( int row_index1, int row_index2, double out_mat[][15])
{
    int n = 15;
    int i;
    double temp[15];

    for(i=0; i<n; i++)
    {
	temp[i] = out_mat[row_index1][i];
	out_mat[row_index1][i] = out_mat[row_index2][i];
	out_mat[row_index2][i] = temp[i];
    }
}

void matexchrow_c16( int row_index1, int row_index2, double out_mat[][16])
{
    int n = 16;
    int i;
    double temp[16];

    for(i=0; i<n; i++)
    {
	temp[i] = out_mat[row_index1][i];
	out_mat[row_index1][i] = out_mat[row_index2][i];
	out_mat[row_index2][i] = temp[i];
    }
}

void matexchrow_c13( int row_index1, int row_index2, double out_mat[][13])
{
    int n = 13;
    int i;
    double temp[13];

    for(i=0; i<n; i++)
    {
	temp[i] = out_mat[row_index1][i];
	out_mat[row_index1][i] = out_mat[row_index2][i];
	out_mat[row_index2][i] = temp[i];
    }
}

void matexchrow_c10( int row_index1, int row_index2, double out_mat[][10])
{
    int n = 10;
    int i;
    double temp[10];

    for(i=0; i<n; i++)
    {
	temp[i] = out_mat[row_index1][i];
	out_mat[row_index1][i] = out_mat[row_index2][i];
	out_mat[row_index2][i] = temp[i];
    }
}

void matexchrow_c6( int row_index1, int row_index2, double out_mat[][6])
{
    int n = 6;
    int i;
    double temp[6];

    for(i=0; i<n; i++)
    {
	temp[i] = out_mat[row_index1][i];
	out_mat[row_index1][i] = out_mat[row_index2][i];
	out_mat[row_index2][i] = temp[i];
    }
}

void matexchrow_c4( int row_index1, int row_index2, double out_mat[][4])
{
    int n = 4;
    int i;
    double temp[4];

    for(i=0; i<n; i++)
    {
	temp[i] = out_mat[row_index1][i];
	out_mat[row_index1][i] = out_mat[row_index2][i];
	out_mat[row_index2][i] = temp[i];
    }
}

void matexchrow_c3( int row_index1, int row_index2, double out_mat[][3])
{
    int n = 3;
    int i;
    double temp[3];

    for(i=0; i<n; i++)
    {
	temp[i] = out_mat[row_index1][i];
	out_mat[row_index1][i] = out_mat[row_index2][i];
	out_mat[row_index2][i] = temp[i];
    }
}

void matsetY( double rotmat[][3], double *c, double out_mat[][6] )
{
    int i, j;

    double lcross[3];
    double skewsym[3][3];

    matmul_c31( 3, rotmat, c, lcross );
    matskewsym( lcross, skewsym );

    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	{
	    if( i==j )
		out_mat[i][j] = 1.;
	    else
		out_mat[i][j] = 0.;

	    out_mat[i][j+3] = -skewsym[i][j];
	}
}

void matsetX( double rotmat[][3], double *r, double out_mat[][6] )
{
    int i, j;

    double lcross[3];
    double skewsym[3][3];

    matmul_c31( 3, rotmat, r, lcross);
    matskewsym( lcross, skewsym );

    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	{
	    if(i==j)
	    {
		out_mat[i][j] = 1.;
		out_mat[i+3][j+3] = 1.;
	    }
	    else
	    {
		out_mat[i][j] = 0.;
		out_mat[i+3][j+3] = 0.;
	    }
	}

    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	    out_mat[i+3][j] = 0.0;

    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	    out_mat[i][j+3] = skewsym[i][j];
}

void matsetinvX( double rotmat[][3], double *r, double out_mat[][6] )
{
    int i, j;

    matsetX( rotmat, r, out_mat );

    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	    out_mat[i][j+3] = -out_mat[i][j+3];
}

void matskewsym( double *vec, double out_mat[][3] )
{
    int i;

    for(i=0; i<3; i++)
	out_mat[i][i] = 0.0;

    out_mat[0][1] = -vec[2];
    out_mat[0][2] =  vec[1];
    out_mat[1][2] = -vec[0];

    out_mat[1][0] = -out_mat[0][1];
    out_mat[2][0] = -out_mat[0][2];
    out_mat[2][1] = -out_mat[1][2];
}

double dabs(double x)
{
	return (x>0?x:-x);
}

void CrossProduct( double *x, double *y, double *out )
{
    out[0] = x[1]*y[2] - y[1]*x[2];
    out[1] = -x[0]*y[2] + x[2]*y[0];
    out[2] = x[0]*y[1] - x[1]*y[0];
}


void RotationMat( int which_axis, double ang, double out_mat[][3])
{
    double cosang = cos(ang);
    double sinang = sin(ang);

    if( which_axis == X_AXIS)
    {
	out_mat[0][0] = 1.;
	out_mat[0][1] = 0.0;
	out_mat[0][2] = 0.0;

	out_mat[1][0] = 0.0;
	out_mat[1][1] = cosang;
	out_mat[1][2] = -sinang;

	out_mat[2][0] = 0.0;
	out_mat[2][1] = sinang;
	out_mat[2][2] = cosang;
    }

    if( which_axis == Y_AXIS )
    {
	out_mat[0][0] = cosang;
	out_mat[0][1] = 0.;
	out_mat[0][2] = sinang;

	out_mat[1][0] = 0.;
	out_mat[1][1] = 1.;
	out_mat[1][2] = 0.;

	out_mat[2][0] = -sinang;
	out_mat[2][1] = 0.;
	out_mat[2][2] = cosang;
    }

    if( which_axis == Z_AXIS )
    {
	out_mat[0][0] = cosang;
	out_mat[0][1] = -sinang;
	out_mat[0][2] = 0.;

	out_mat[1][0] = sinang;
	out_mat[1][1] = cosang;
	out_mat[1][2] = 0.;

	out_mat[2][0] = 0.;
	out_mat[2][1] = 0.;
	out_mat[2][2] = 1.;
    }
}

void  RotationMatT( int which_axis, double ang, double out_mat[][4])
{
    double cosang = cos(ang);
    double sinang = sin(ang);

    if( which_axis == X_AXIS)
    {
	out_mat[0][0] = 1.;
	out_mat[0][1] = 0.0;
	out_mat[0][2] = 0.0;

	out_mat[1][0] = 0.0;
	out_mat[1][1] = cosang;
	out_mat[1][2] = -sinang;

	out_mat[2][0] = 0.0;
	out_mat[2][1] = sinang;
	out_mat[2][2] = cosang;
    }

    if( which_axis == Y_AXIS )
    {
	out_mat[0][0] = cosang;
	out_mat[0][1] = 0.;
	out_mat[0][2] = sinang;

	out_mat[1][0] = 0.;
	out_mat[1][1] = 1.;
	out_mat[1][2] = 0.;

	out_mat[2][0] = -sinang;
	out_mat[2][1] = 0.;
	out_mat[2][2] = cosang;
    }

    if( which_axis == Z_AXIS )
    {
	out_mat[0][0] = cosang;
	out_mat[0][1] = -sinang;
	out_mat[0][2] = 0.;

	out_mat[1][0] = sinang;
	out_mat[1][1] = cosang;
	out_mat[1][2] = 0.;

	out_mat[2][0] = 0.;
	out_mat[2][1] = 0.;
	out_mat[2][2] = 1.;
    }


	out_mat[0][3] = 0.0;
	out_mat[1][3] = 0.0;
	out_mat[2][3] = 0.0;

    out_mat[3][0] = 0.0;
    out_mat[3][1] = 0.0;
    out_mat[3][2] = 0.0;
    out_mat[3][3] = 1.0;

}

void Get3OrdPolyCoeff( double t_0, double t_f, double y_0, double y_f, double *out_coeff )
{
    double mat[4][4];
    double vec[4];
    double invmat[4][4];

    mat[0][0] = 1.0;
    mat[0][1] = t_0;
    mat[0][2] = t_0*t_0;
    mat[0][3] = t_0*t_0*t_0;

    mat[1][0] = 1.;
    mat[1][1] = t_f;
    mat[1][2] = t_f*t_f;
    mat[1][3] = t_f*t_f*t_f;

    mat[2][0] = 0.;
    mat[2][1] = 1.;
    mat[2][2] = 2.*t_0;
    mat[2][3] = 3.*t_0*t_0;

    mat[3][0] = 0.0;
    mat[3][1] = 1.;
    mat[3][2] = 2.*t_f;
    mat[3][3] = 3.*t_f*t_f;

    vec[0] = y_0;
    vec[1] = y_f;
    vec[2] = 0.;
    vec[3] = 0.;

    matinv_c4( mat, invmat );
    matmul_c41( 4, invmat, vec, out_coeff);
}

void  Get3OrdPolyCoeff_b( double t_f, double y_0, double ydot_0, double y_f, double *out_coeff )
{
    /*
    double mat[2][2];
    double vec[2], vec2[2];
    double invmat[2][2];

    out_coeff[0] = y_0;
    out_coeff[1] = ydot_0;

    vec[0] = -out_coeff[0]-out_coeff[1]*t_f + y_f;
    vec[1] = -out_coeff[1];

    mat[0][0] = t_f*t_f;
    mat[0][1] = t_f*t_f*t_f;
    mat[1][0] = 2.*t_f;
    mat[1][1] = 3.*t_f*t_f;

    matinv_c2( mat, invmat );
    matmul_c21( 2, invmat, vec, vec2 );
    out_coeff[2] = vec2[0];
    out_coeff[3] = vec2[1];
    */

    out_coeff[0] = y_0;
    out_coeff[1] = ydot_0;
    out_coeff[2] = (3.*(-out_coeff[0]-out_coeff[1]*t_f+y_f)+out_coeff[1]*t_f)/t_f/t_f;
    out_coeff[3] = (-out_coeff[1]*t_f-2.*(-out_coeff[0]-out_coeff[1]*t_f+y_f))/t_f/t_f/t_f;
}

void Get2OrdRoot( double coeff_a, double coeff_b, double coeff_c, double *out_root1, double *out_root2)
{
    double temp1, temp2;

    temp1 = -coeff_b/2./coeff_a;
    temp2 = sqrt((coeff_b/2./coeff_a)*(coeff_b/2./coeff_a) - coeff_c/coeff_a);

    *out_root1 = temp1 + temp2;
    *out_root2 = temp1 - temp2;
}

double Tsinh( double x )
{
    double ret;
    double yi = x;
    int i;

    ret = yi;
    for(i=2; i<= TAYLOR_NUM_SINH; i++)
    {
	yi *= x*x/(double)((2*i-1)*(2*i-2));
	ret += yi;
    }

    return ret;
}

double Tcosh( double x )
{
    double ret;
    double yi = 1.;
    int i;

    ret = 1.;
    for(i=2; i<= TAYLOR_NUM_COSH; i++)
    {
	yi *= x*x/(double)((2*i-2)*(2*i-3));
	ret += yi;
    }

    return ret;
}

double Tln( double x )
{
    int n = 0;
    double A = x;
    double ret;
    double yi;
    int i;

    while( A > 2. )
    {
        A /= (double)NATURAL_E;
	n++;
    }

    yi = A-1.;
    ret = yi;
    for(i=2; i<=TAYLOR_NUM_LN; i++)
    {
	yi *= (double)((-1)*(i-1))/(double)(i)*(A-1.0);
	ret += yi;
    }

    ret += (double)n;

    return ret;
}

void  GetQuaternion( double Rot[][3], double out_Q[4] )
{
    double temp[4];
    int i;

    temp[0] = Rot[0][0] + Rot[1][1] + Rot[2][2] + 1.;
    temp[1] = Rot[0][0] - Rot[1][1] - Rot[2][2] + 1.;
    temp[2] = Rot[1][1] - Rot[2][2] - Rot[0][0] + 1.;
    temp[3] = Rot[2][2] - Rot[0][0] - Rot[1][1] + 1.;

    for(i=0; i<4; i++)
    {
	if(temp[i] <= 0.0)
	    temp[i] = 0.0;
    }

    out_Q[0] = 1./2. * sqrt(temp[0]);
    out_Q[1] = 1./2. * mysgn(Rot[2][1] - Rot[1][2]) * sqrt(temp[1]);
    out_Q[2] = 1./2. * mysgn(Rot[0][2] - Rot[2][0]) * sqrt(temp[2]);
    out_Q[3] = 1./2. * mysgn(Rot[1][0] - Rot[0][1]) * sqrt(temp[3]);
}

double mysgn( double value )
{
    if ( value > 0. )
	return 1.;
    else
	return -1.;
}

double Qt_Norm( double *qt )
{
    double ret;
    int i;

    ret = 0.0;
    for(i=0; i<4; i++)
	ret += qt[i]*qt[i];

    return ret;
}

void Qt_ToMatrix( double *qt, double out_M[][3] )
{
    double Nq = Qt_Norm( qt );
    double s = (Nq > 0.0) ? (2.0/Nq) : 0.0 ;
    double xs = qt[1]*s,   ys = qt[2]*s,   zs = qt[3]*s;
    double wx = qt[0]*xs,  wy = qt[0]*ys,  wz = qt[0]*zs;
    double xx = qt[1]*xs,  xy = qt[1]*ys,  xz = qt[1]*zs;
    double yy = qt[2]*ys,  yz = qt[2]*zs,  zz = qt[3]*zs;
    out_M[0][0] = 1.0 - (yy+zz);
    out_M[1][0] = xy + wz;
    out_M[2][0] = xz - wy;
    out_M[0][1] = xy - wz;
    out_M[1][1] = 1.0 - (xx+zz);
    out_M[2][1] = yz + wx;
    out_M[0][2] = xz + wy;
    out_M[1][2] = yz - wx;
    out_M[2][2] = 1.0 - (xx+yy);
}

void  Qt_FromMatrix( double M[][3], double *out_q )
{
    double tr, s;

    tr = M[0][0] + M[1][1] + M[2][2];
    if (tr >= 0.0)
    {
	s = sqrt(tr + 1.0);
	out_q[0] = s*0.5;
	s = 0.5/s;
	out_q[1] = (M[2][1] - M[1][2])*s;
	out_q[2] = (M[0][2] - M[2][0])*s;
	out_q[3] = (M[1][0] - M[0][1])*s;
    }
    else
    {
	int h = 0;
	if( M[1][1] > M[0][0] )
	    h = 1;
	if( M[2][2] > M[h][h] )
	    h = 2;
	switch(h)
	{
	    case 0 :
		caseMacro(0, 1, 2, 0, 1, 2, M, out_q);
		break;
	    case 1:
		caseMacro(1, 2, 0, 1, 2, 0, M, out_q);
		break;
	    case 2:
		caseMacro(2, 0, 1, 2, 0, 1, M, out_q);
		break;
	}
    }
}

void  caseMacro(int i, int j, int k, int I, int J, int K, double M[][3], double *out_q)
{
    double s;

    s = sqrt( (M[I][I] - (M[J][J]+M[K][K])) + 1.0 );
    out_q[i+1] = s*0.5;
    s = 0.5/s;
    out_q[j+1] = (M[I][J] + M[J][I])*s;
    out_q[k+1] = (M[K][I] + M[I][K])*s;
    out_q[0]   = (M[K][J] - M[J][K])*s;
}

void  rpy_FromMatrix( double R[][3], double *out_v )
{
    double alpha, beta, gamma;

    alpha = atan2( R[1][0], R[0][0] );
    beta = atan2( -R[2][0], sqrt(R[0][0]*R[0][0] + R[1][0]*R[1][0]) );
    //beta = atan2( -R[2][0], sqrt(R[2][1]*R[2][1] + R[2][2]*R[2][2]) );
    gamma = atan2( R[2][1], R[2][2] );

    out_v[0] = gamma;	    // roll
    out_v[1] = beta;	    // pitch
    out_v[2] = alpha;	    // yaw
}

void  AxisAngle_FromMatrix( double R[][3], double *k, double *theta )
{
	double angle;

	angle = ArcCos( ( R[0][0] + R[1][1] + R[2][2] -1.0 ) / 2.0 );
	if( angle != 0.0 ){
		k[0] = ( R[2][1] - R[1][2] )/( 2.0*sin(angle) );
		k[1] = ( R[0][2] - R[2][0] )/( 2.0*sin(angle) );
		k[2] = ( R[1][0] - R[0][1] )/( 2.0*sin(angle) );
	} else {
		k[0] = k[1] = k[2] = 1.0/sqrt(3.0);
	}
	*theta = angle;
}


void  Euler_XYZ_FromMatrix( double R[][3], double *out_v )
{
    double alpha, beta, gamma;

    alpha = atan2( -R[1][2], R[2][2] );
    beta  = atan2(  R[0][2], sqrt(R[0][0]*R[0][0] + R[0][1]*R[0][1]) );
    gamma = atan2( -R[0][1], R[0][0] );

    out_v[0] = alpha;
    out_v[1] = beta;
    out_v[2] = gamma;
}

void  Get_R_from_AngleAxis( double *k, double theta, double R[][3] )
{
	double kx, ky, kz;
	double vt, ct, st;
	kx = k[0];
	ky = k[1];
	kz = k[2];

	vt = 1-cos(theta);
	ct = cos(theta);
	st = sin(theta );

	R[0][0] = kx*kx*vt + ct;
	R[0][1] = kx*ky*vt - kz*st;
	R[0][2] = kx*kz*vt + ky*st;
	R[1][0] = kx*ky*vt + kz*st;
	R[1][1] = ky*ky*vt + ct;
	R[1][2] = ky*kz*vt - kx*st;
	R[2][0] = kx*kz*vt - ky*st;
	R[2][1] = ky*kz*vt + kx*st;
	R[2][2] = kz*kz*vt + ct;
}

void Get_R_from_rpy( double *rpy_vec, double out_R[][3] )
{
	double gamma, beta, alpha;
    double sinalpha, cosalpha, sinbeta, cosbeta, singamma, cosgamma;

	gamma = rpy_vec[0]; // roll ( X )
	beta  = rpy_vec[1]; // pitch ( Y )
	alpha = rpy_vec[2]; // yaw ( Z )

    sinalpha = sin(alpha);
    cosalpha = cos(alpha);
    sinbeta  = sin(beta );
    cosbeta  = cos(beta );
	singamma = sin(gamma);
	cosgamma = cos(gamma);

	out_R[0][0] =  cosalpha*cosbeta;
	out_R[0][1] = -sinalpha*cosgamma + cosalpha*sinbeta*singamma;
    out_R[0][2] =  sinalpha*singamma + cosalpha*sinbeta*cosgamma;
    out_R[1][0] =  sinalpha*cosbeta;
    out_R[1][1] =  cosalpha*cosgamma + sinalpha*sinbeta*singamma;
	out_R[1][2] = -cosalpha*singamma + sinalpha*sinbeta*cosgamma;
    out_R[2][0] = -sinbeta;
    out_R[2][1] =  cosbeta*singamma;
    out_R[2][2] =  cosbeta*cosgamma;
}

void Get_T_from_rpy( double *rpy_vec, double out_T[][4] )
{
	double gamma, beta, alpha;
    double sinalpha, cosalpha, sinbeta, cosbeta, singamma, cosgamma;

	gamma = rpy_vec[0]; // roll ( X )
	beta  = rpy_vec[1]; // pitch ( Y )
	alpha = rpy_vec[2]; // yaw ( Z )

    sinalpha = sin(alpha);
    cosalpha = cos(alpha);
    sinbeta  = sin(beta );
    cosbeta  = cos(beta );
	singamma = sin(gamma);
	cosgamma = cos(gamma);

	out_T[0][0] =  cosalpha*cosbeta;
	out_T[0][1] = -sinalpha*cosgamma + cosalpha*sinbeta*singamma;
    out_T[0][2] =  sinalpha*singamma + cosalpha*sinbeta*cosgamma;
    out_T[1][0] =  sinalpha*cosbeta;
    out_T[1][1] =  cosalpha*cosgamma + sinalpha*sinbeta*singamma;
	out_T[1][2] = -cosalpha*singamma + sinalpha*sinbeta*cosgamma;
    out_T[2][0] = -sinbeta;
    out_T[2][1] =  cosbeta*singamma;
    out_T[2][2] =  cosbeta*cosgamma;
	out_T[3][0] = 0.0;
	out_T[3][1] = 0.0;
	out_T[3][2] = 0.0;

	out_T[0][3] = 0.0;
	out_T[1][3] = 0.0;
	out_T[2][3] = 0.0;
	out_T[3][3] = 1.0;
}

void  Get_T_from_Euler_ZYX( double *vec, double out_T[][3] )
{
    double sinalpha, cosalpha, sinbeta, cosbeta, singamma, cosgamma;

    sinalpha = sin(vec[0]);
    cosalpha = cos(vec[0]);
    sinbeta  = sin(vec[1]);
    cosbeta  = cos(vec[1]);
	singamma = sin(vec[2]);
	cosgamma = cos(vec[2]);

	out_T[0][0] =  cosalpha*cosbeta;
	out_T[0][1] = -sinalpha*cosgamma + cosalpha*sinbeta*singamma;
    out_T[0][2] =  sinalpha*singamma + cosalpha*sinbeta*cosgamma;
    out_T[1][0] =  sinalpha*cosbeta;
    out_T[1][1] =  cosalpha*cosgamma + sinalpha*sinbeta*singamma;
	out_T[1][2] = -cosalpha*singamma + sinalpha*sinbeta*cosgamma;
    out_T[2][0] = -sinalpha;
    out_T[2][1] =  cosbeta*singamma;
    out_T[2][2] =  cosbeta*cosgamma;
}


void  Get_T_from_Euler_XYZ( double *vec, double out_T[][3] )
{
    double sinalpha, cosalpha, sinbeta, cosbeta;

    sinalpha = sin(vec[0]);
    cosalpha = cos(vec[0]);
    sinbeta  = sin(vec[1]);
    cosbeta  = cos(vec[1]);

    out_T[0][0] = 1.0;
    out_T[0][1] = 0.0;
    out_T[0][2] = sinbeta;
    out_T[1][0] = 0.0;
    out_T[1][1] = cosalpha;
    out_T[1][2] = -sinalpha*cosbeta;
    out_T[2][0] = 0.0;
    out_T[2][1] = sinalpha;
    out_T[2][2] = cosalpha*cosbeta;
}

void  Get_T_from_Rot( double R[][3], double *out_Euler_XYZ, double out_T[][3] )
{
    double alpha, beta, gamma;
    double sinalpha, cosalpha, sinbeta, cosbeta;

    alpha = atan2( R[1][2], -R[2][2] );
    beta  = atan2( R[0][2], -sqrt(R[0][0]*R[0][0] + R[0][1]*R[0][1]) );
    gamma = atan2( R[0][1], -R[0][0] );

    sinalpha = sin(alpha);
    sinbeta  = sin(beta);
    cosalpha = cos(alpha);
    cosbeta  = cos(beta);

    out_Euler_XYZ[0] = alpha;
    out_Euler_XYZ[1] = beta;
    out_Euler_XYZ[2] = gamma;

    out_T[0][0] = 1.;
    out_T[0][1] = 0.;
    out_T[0][2] = sinbeta;
    out_T[1][0] = 0.;
    out_T[1][1] = cosalpha;
    out_T[1][2] = -sinalpha*cosbeta;
    out_T[2][0] = 0.;
    out_T[2][1] = sinalpha;
    out_T[2][2] = cosalpha*cosbeta;
}

double in_prod( int dim, double *v1, double *v2 )
{
    double ret;
    int i;

    ret = 0.;

    for(i=0; i<dim; i++)
	ret += v1[i]*v2[i];

    return ret;
}

double vsum( int dim, double *v )
{
	int i;
    double sum;
	sum=0;
	for(i=0; i<dim; i++ ){
		sum += v[i];
	}

    return sum;
}


double vnorm( int dim, double *v )
{
    double inprod;

    inprod = in_prod( dim, v, v );

    return sqrt(inprod);
}

void vnormalize( int dim, double *v )
{
    double inprod;
	int i;

    inprod = in_prod( dim, v, v );
    inprod = sqrt(inprod);

	for(i=0; i<dim; i++) v[i] /= inprod;
}


double ArcCos( double val )	// output : 0 ~ PI
{
    double costheta, sintheta;

    costheta = val;
    if(costheta > 1.)
		costheta = 1.0;
    if(costheta < -1.)
		costheta = -1.0;

    sintheta = sqrt( 1. - costheta*costheta );

    return atan2( sintheta, costheta );
}

double ArcSin( double val )	// output : -PI/2 ~ PI/2
{
    double costheta, sintheta;

    sintheta = val;
    if(sintheta > 1.)
	sintheta = 1.0;
    if(sintheta < -1.)
	sintheta = -1.0;
    costheta = sqrt( 1. - sintheta*sintheta );

    return atan2( sintheta, costheta );
}

double ArcTan( double length_bottom, double length_height )
{
    double hypotenuse;
    double costheta, sintheta;

    hypotenuse = sqrt(length_bottom*length_bottom + length_height*length_height);

    costheta = length_bottom/hypotenuse;
    sintheta = length_height/hypotenuse;

    return atan2( sintheta, costheta );
}

void  Get_mD( double mass, double c_part[3], double c_tot[3], double out_I[][3] )	// kg m^2
{
    int i, j;
    double cvec[3], kcvec[3];
    double skew[3][3], kskew[3][3];

    matsub_c1( 3, c_part, c_tot, cvec );
    for(i=0; i<3; i++)
	kcvec[i] = cvec[i];

    matskewsym( cvec, skew );
    matskewsym( kcvec, kskew );

    matmul_c33( 3, skew, kskew, out_I);

    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	    out_I[i][j] *= -mass*out_I[i][j];
}

void  Get_RotI( double Inertia[][3], double Rot[][3], double out_I[][3] )
{
    int i, j;
    double Transpose_of_Rot[3][3];

    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	    Transpose_of_Rot[i][j] = Rot[j][i];

    matmul_c33( 3, Rot, Inertia, out_I );
    matrmul_r33c33( Transpose_of_Rot, out_I);
}

void  Set_Homo_Trans_from_Rot( double Rot[][3], double out_Homo_Trans[][4] )
{
    int i, j;

    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
		out_Homo_Trans[i][j] = Rot[i][j];
}

void  Get_Rot_from_Homo_Trans( double Homo_Trans[][4], double out_Rot[][3] )
{
    int i, j;

    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	    out_Rot[i][j] = Homo_Trans[i][j];
}


double GetDist_2D( double *vecA, double *vecB )
{
    int i;
    int dim = 2;
    double vecAB[2];

    for(i=0; i<dim; i++)
	vecAB[i] = vecB - vecA;

    return GetNorm_2D( vecAB );
}

double GetNorm_2D( double *vec )
{
    int dim = 2;

    int i;
    double sum = 0.0;

    for(i=0; i<dim; i++)
	sum += vec[i]*vec[i];

    return (double)sqrt(sum);
}

void  GetUnitDirVec_2D( double *outvec, double *vecStart, double *vecEnd )
{
    int dim = 2;
    int i;
    double Result[2];
    double norm;

    for(i=0; i<dim; i++)
	Result[i] = vecEnd[i] - vecStart[i];

    norm = GetNorm_2D( Result );

    for(i=0; i<dim; i++)
	Result[i] = Result[i]/norm;

    for(i=0; i<dim; i++)
	outvec[i] = Result[i];
}

void  SetRotationMat( double ang, double out_mat[][2])
{
    double cosang = cos(ang);
    double sinang = sin(ang);

    out_mat[0][0] = cosang;
    out_mat[0][1] = -sinang;
    out_mat[1][0] = sinang;
    out_mat[1][1] = cosang;
}


void T_from_AxisAnglePosition( double *in, double T[][4] )
{
	double R[3][3];
	int i;

	Get_R_from_AngleAxis( &(in[3]), in[6] , R);
	matI_c4(T);
	Set_Homo_Trans_from_Rot( R, T );
	for(i=0; i<3; i++)
		T[i][3] = in[i];

}

void PositionAxisAngle_from_T( double T[][4], double *out )
{
	double R[3][3];
	double k[3], theta;
	int i;

	Get_Rot_from_Homo_Trans(T, R );
	AxisAngle_FromMatrix(R, k, &theta );

	for( i=0; i<3; i++){
		out[i] = T[i][3];
		out[i+3] = k[i];
	}
	out[6] = theta;
}
