#ifndef _MATH_LIB_H_
#define _MATH_LIB_H_


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

	#ifndef M_PI_2
	#define M_PI_2 1.57079632679489661922
	#endif

	#ifndef M_PI
	#define M_PI 3.14159265358979323846
	#endif

	#define X_AXIS   1
	#define Y_AXIS   2
	#define Z_AXIS   3

	#define TAYLOR_NUM_SINH	    15
	#define TAYLOR_NUM_COSH	    15
	#define TAYLOR_NUM_LN	    40

	#define NATURAL_E	    2.7182818284590
	#define SINGULAR_TOLERANCE 0.001

	void matsavefile(char *str, double **M, int m, int n);
	void matsavefilex(char *str, double **M, int ms, int mf, int ns, int nf);

/***************************************************************
 * matrix functions for general case
****************************************************************/
	// matrix & vector allocation
	double **smatrix(unsigned int m, unsigned int n);
	double *svector(unsigned int m);

	double **dmatrix(int nrl, int nrh, int ncl, int nch);
	double *dvector(int nl, int nh);

	// SVD functions svd.c
	void svdcmp(double **a, int m, int n, double w[], double **v);

	// matrix svd pseudo inversion with fixed sinuglar tollerance
	unsigned int matsvdinvall( double **mat, int m, int n, double singular_tollerance, double *det, double **out_mat);
	unsigned int matsvdinv( double **mat, int m, int n, double **out_mat);
	unsigned int matsvdinvDLS( double **mat, int m, int n, double lamda, double **out_mat);

	// trnaspose
	void mattranspose( double **mat, int m, int n, double **out_mat);

	// multiplication
	void matmul_vec(double **mat, int m, int n, double *rvec, double *out_vec );
	void matmul_cc( int lm, int ln, int rn, double **lmat, double **rmat, double **out_mat);

	// null space of matrix
	void matnull( int m, int n, double **J, double **Ji, double **N);


/***************************************************************
 * matrix functions for fixed dimension
****************************************************************/
	void matzero_c4( int m, double out_mat[][4] );
	void matzero_c6( int m, double out_mat[][6] );
	void matzero_c1( int m, double *out_vec );
	void matzero_c7( int m, double out_mat[][7] );

	void matmul_c42( int lm, double lmat[][4], double rmat[][2], double out_mat[][2]);
	void matmul_c44( int lm, double lmat[][4], double rmat[][4], double out_mat[][4]);
	void matmul_c36( int lm, double lmat[][3], double rmat[][6], double out_mat[][6]);
	void matmul_c63( int lm, double lmat[][6], double rmat[][3], double out_mat[][3]);
	void matmul_c613( int lm, double lmat[][6], double rmat[][13], double out_mat[][13]);
	void matmul_c1313( int lm, double lmat[][13], double rmat[][13], double out_mat[][13]);
	void matmul_c1515( int lm, double lmat[][15], double rmat[][15], double out_mat[][15]);
	void matmul_c1616( int lm, double lmat[][16], double rmat[][16], double out_mat[][16]);
	void matmul_c133( int lm, double lmat[][13], double rmat[][3], double out_mat[][3]);
	void matmul_c153( int lm, double lmat[][15], double rmat[][3], double out_mat[][3]);
	void matmul_c164( int lm, double lmat[][16], double rmat[][4], double out_mat[][4]);
	void matmul_c313( int lm, double lmat[][3], double rmat[][13], double out_mat[][13]);
	void matmul_c316( int lm, double lmat[][3], double rmat[][16], double out_mat[][16]);
	void matmul_c37( int lm, double lmat[][3], double rmat[][7], double out_mat[][7]);
	void matmul_c34( int lm, double lmat[][3], double rmat[][4], double out_mat[][4]);
	void matmul_c66( int lm, double lmat[][6], double rmat[][6], double out_mat[][6]);
	void matmul_c33( int lm, double lmat[][3], double rmat[][3], double out_mat[][3]);

	void matmul_c21( int lm, double lmat[][2], double *rvec, double *out_vec);
	void matmul_c31( int lm, double lmat[][3], double *rvec, double *out_vec);
	void matmul_c61( int lm, double lmat[][6], double *rvec, double *out_vec);
	void matmul_c71( int lm, double lmat[][7], double *rvec, double *out_vec);
	void matmul_c101( int lm, double lmat[][10], double *rvec, double *out_vec);
	void matmul_c131( int lm, double lmat[][13], double *rvec, double *out_vec);
	void matmul_c161( int lm, double lmat[][16], double *rvec, double *out_vec);
	void matmul_c41( int lm, double lmat[][4], double *rvec, double *out_vec);
	void matmul_c( int row, int col,  double **lmat, double *rvec, double *out_vec);


	// rightside matrix multiplication
	void matrmul_r66c66( double rmat[][6], double out_mat[][6]);
	void matrmul_r44c44( double rmat[][4], double out_mat[][4]);
	void matrmul_r33c33( double rmat[][3], double out_mat[][3]);

	// leftside matrix multiplication
	void matlmul_r44c44( double lmat[][4], double out_mat[][4]);
	void matlmul_r33c36( double lmat[][3], double out_mat[][6]);

	// scalar multiplication
	void matsmul_c7( int m, double scalar, double out_mat[][7]);
	void matsmul_c6( int m, double scalar, double out_mat[][6]);
	void matsmul_c4( int m, double scalar, double out_mat[][4]);
	void matsmul_c1( int m, double scalar, double *out_vec);

	// addition
	void matadd_c7( int m, double mat[][7], double out_mat[][7]);
	void matadd_c4( int m, double mat[][4], double out_mat[][4]);
	void matadd_c3( int m, double mat[][3], double out_mat[][3]);
	void matadd_c6( int m, double mat[][6], double out_mat[][6]);
	void matadd_c1( int m, double *vec, double *out_vec);
	void matadd_c( int m, double *lvec, double *rvec, double *out_vec );

	// subtraction
	void matsub_c1( int m, const double *lvec, const double *rvec, double *out_vec );

	// transpose
	void mattransp_r3c6( double mat[][6], double out_mat[][3] );
	void mattransp_r3c13( double mat[][13], double out_mat[][3] );
	void mattransp_r3c15( double mat[][15], double out_mat[][3] );
	void mattransp_r4c16( double mat[][16], double out_mat[][4] );

	// set identity matrix
	void matI_c13( double out_mat[][13]);
	void matI_c16( double out_mat[][16]);
	void matI_c15( double out_mat[][15]);
	void matI_c10( double out_mat[][10]);
	void matI_c6( double out_mat[][6]);
	void matI_c4( double out_mat[][4]);
	void matI_c3( double out_mat[][3]);

	// copy matrix
	void matcp_c10( int m, double mat[][10], double out_mat[][10] );
	void matcp_c13( int m, double mat[][13], double out_mat[][13] );
	void matcp_c16( int m, double mat[][16], double out_mat[][16] );
	void matcp_c15( int m, double mat[][15], double out_mat[][15] );
	void matcp_c7( int m, double mat[][7], double out_mat[][7] );
	void matcp_c6( int m, double mat[][6], double out_mat[][6] );
	void matcp_c4( int m, double mat[][4], double out_mat[][4] );
	void matcp_c3( int m, double mat[][3], double out_mat[][3] );
	void matcp_c1( int m, double *vec, double *out_vec);

	// matrix inversion
	void matinv_c2( double mat[][2], double out_mat[][2] );
	void matinv_c3( double mat[][3], double out_mat[][3] );
	void matinv_c4( double mat[][4], double out_mat[][4] );
	void matinv_c6( double mat[][6], double out_mat[][6] );
	void matinv_c10( double mat[][10], double out_mat[][10] );
	void matinv_c13( double mat[][13], double out_mat[][13] );
	void matinv_c15( double mat[][15], double out_mat[][15] );
	void matinv_c16( double mat[][16], double out_mat[][16] );

	// matrix pseudo inversion
	void matpinv_r3c6( double mat[][6], double weight[][6], double out_mat[][3]);
	void matpinv_r3c13( double mat[][13], double weight[][13], double out_mat[][3]);
	void matpinv_r3c15( double mat[][15], double weight[][15], double out_mat[][3]);
	void matpinv_r4c16( double mat[][16], double weight[][16], double out_mat[][4]);
	void matpinv_r1c13( double* row_vec, double weight[][13], double* out_col_vec);

	// private partial pivoting method
	int matpivot_c3( int nm, double mat[][3] );
	int matpivot_c4( int nm, double mat[][4] );
	int matpivot_c6( int nm, double mat[][6] );
	int matpivot_c10( int nm, double mat[][10] );
	int matpivot_c13( int nm, double mat[][13] );
	int matpivot_c16( int nm, double mat[][16] );
	int matpivot_c15( int nm, double mat[][15] );

	// exchange two rows in a matrix
	void matexchrow_c3( int row_index1, int row_index2, double out_mat[][3]);
	void matexchrow_c4( int row_index1, int row_index2, double out_mat[][4]);
	void matexchrow_c6( int row_index1, int row_index2, double out_mat[][6]);
	void matexchrow_c10( int row_index1, int row_index2, double out_mat[][10]);
	void matexchrow_c13( int row_index1, int row_index2, double out_mat[][13]);
	void matexchrow_c16( int row_index1, int row_index2, double out_mat[][16]);
	void matexchrow_c15( int row_index1, int row_index2, double out_mat[][15]);

	// skew symmetric matrix (cross product)
	void matskewsym( double *vec, double out_mat[][3] );

	void matsetX( double rotmat[][3], double *r, double out_mat[][6] );
	void matsetinvX( double rotmat[][3], double *r, double out_mat[][6] );

	void matsetY( double rotmat[][3], double *c, double out_mat[][6] );



/***************************************************************
 * matrix utils
****************************************************************/
	double dabs(double x);
	void  CrossProduct( double *x, double *y, double *out);
	void  RotationMat( int which_axis, double ang, double out_mat[][3]);
	void  RotationMatT( int which_axis, double ang, double out_mat[][4]);
	void  SetRotationMat( double ang, double out_mat[][2]);
	void  Get3OrdPolyCoeff( double t_0, double t_f, double y_0, double y_f, double *out_coeff );
	void  Get3OrdPolyCoeff_b( double t_f, double y_0, double ydot_0, double y_f, double *out_coeff );
	void  Get2OrdRoot( double coeff_a, double coeff_b, double coeff_c, double *out_root1, double *out_root2);
	double Tsinh( double x );
	double Tcosh( double x );
	double Tln( double x );
	void  GetQuaternion( double Rot[][3], double out_Q[4] );
	double mysgn( double value );
	double Qt_Norm( double *qt );
	void  Qt_ToMatrix( double *qt, double out_M[][3] );
	void  Qt_FromMatrix( double M[][3], double *out_q );
	void  caseMacro(int i, int j, int k, int I, int J, int K, double M[][3], double *out_q);
	void  rpy_FromMatrix( double R[][3], double *out_v );
	void  AxisAngle_FromMatrix( double R[][3], double *k, double *theta );
	void  Euler_XYZ_FromMatrix( double R[][3], double *out_v );
	void  Get_R_from_AngleAxis( double *k, double theta, double R[][3] );
	void  Get_R_from_rpy( double *rpy_vec, double out_R[][3] );
	void  Get_T_from_rpy( double *rpy_vec, double out_T[][4] );
	void  Get_T_from_Euler_ZYX( double *vec, double out_T[][3] );
	void  Get_T_from_Euler_XYZ( double *vec, double out_T[][3] );
	void  Get_T_from_Rot( double R[][3], double *out_Euler_XYZ, double out_T[][3] );
	double in_prod( int dim, double *v1, double *v2 );
	double ArcCos( double val );
	double ArcSin( double val );
	double ArcTan( double length_bottom, double length_height );

	double vsum( int dim, double *v );
	double vnorm( int dim, double *v );
	void   vnormalize( int dim, double *v );

	void  Get_mD( double mass, double c_part[3], double c_tot[3], double out_I[][3] );
	void  Get_RotI( double Inertia[][3], double Rot[][3], double out_I[][3] );
	void  Set_Homo_Trans_from_Rot( double Rot[][3], double out_Homo_Trans[][4] );
	void  Get_Rot_from_Homo_Trans( double Homo_Trans[][4], double out_Rot[][3] );
	double GetDist_2D( double *vecA, double *vecB );
	double GetNorm_2D( double *vec );
	void  GetUnitDirVec_2D( double *outvec, double *vecStart, double *vecEnd );
	void T_from_AxisAnglePosition( double *in, double T[][4] );
	void PositionAxisAngle_from_T( double T[][4], double *out );
	//void Homo_Inverse( double Homo_in[][4], double Homo_out[][4]);


#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif //_MATH_LIB_H_
