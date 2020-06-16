/*******************************************************************************
    This file contains the cylindrification routines from HARMPI,
    originally written by Sasha Tchekhovskoy.

    It also contains routines to transform coordinates from HARMPI's
    cylindified MKS (here called CMKS) into standard MKS and then
    Boyer-Linquist coordinates used in grmonty.

    Some of these routines are based on Jason Dexter's grtrans code,
    others are based on EHT-babel.
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "definitions.h"
#include "declarations.h"


void read_metric(char *fname)
{
    int i, j, k;
    int l, m, n;
    int index;
    double var[96];
    FILE *fp;

    fp = fopen(fname, "rb");

    if(fp==NULL) {
      fprintf(stderr,"error opening metricfile\n") ;
      exit(1) ;
    }

    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {

                fread(var, sizeof(double), 96, fp);
                index = 0;

                for (m = 0; m < NDIM; m++) {
                    for (n = 0; n < NDIM; n++) {
                        grid_gcov[i][j][k][m][n] = var[index];
                        index++;
                    }
                }

                for (m = 0; m < NDIM; m++) {
                    for (n = 0; n < NDIM; n++) {
                        grid_gcon[i][j][k][m][n] = var[index];
                        index++;
                    }
                }

                for (l = 0; l < NDIM; l++) {
                    for (m = 0; m < NDIM; m++) {
                        for (n = 0; n < NDIM; n++) {
                            grid_conn[i][j][k][l][m][n] = var[index];
                            index++;
                        }
                    }
                }

            }
        }
    }
    fclose(fp);
}


void write_metric(char *fname)
{
    int i, j, k;
    int l, m, n;
    FILE *fp;

    fp = fopen(fname, "wb");

    if(fp==NULL) {
      fprintf(stderr,"error opening metricfile\n") ;
      exit(1) ;
    }

    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {

                for (m = 0; m < NDIM; m++) {
                    for (n = 0; n < NDIM; n++) {
                        fwrite(&grid_gcov[i][j][k][m][n], sizeof(double), 1, fp);
                    }
                }

                for (m = 0; m < NDIM; m++) {
                    for (n = 0; n < NDIM; n++) {
                        fwrite(&grid_gcon[i][j][k][m][n], sizeof(double), 1, fp);
                    }
                }

                for (l = 0; l < NDIM; l++) {
                    for (m = 0; m < NDIM; m++) {
                        for (n = 0; n < NDIM; n++) {
                            fwrite(&grid_conn[i][j][k][l][m][n], sizeof(double), 1, fp);
                        }
                    }
                }

            }
        }
    }
    fclose(fp);
}


// insert metric here
void gcov_func(double *X, double gcovp[][NDIM])
{
  int i,j,k,l ;
  double sth,cth,s2,rho2 ;
  double rr,tth,pphi ;
  double tfac,rfac,hfac,pfac ;
  double ggcov[NDIM][NDIM];
  double ddxdxp[NDIM][NDIM];

  DLOOP ggcov[j][k] = 0. ;
  DLOOP gcovp[j][k] = 0.;

  bl_coord(X,&rr,&tth,&pphi) ;

  cth = cos(tth) ;
  sth = sin(tth) ;

  s2 = sth*sth ;
  rho2 = rr*rr + a*a*cth*cth ;

  //compute Jacobian x1,x2,x3 -> r,th,phi
  dxdxp_func(X, ddxdxp);

  if(GR==1){
    ggcov[0][0] = (-1. + 2.*rr/rho2);
    ggcov[0][1] = (2.*rr/rho2);
    ggcov[0][3] = (-2.*a*rr*s2/rho2);

    ggcov[1][0] = ggcov[0][1] ;
    ggcov[1][1] = (1. + 2.*rr/rho2);
    ggcov[1][3] = (-a*s2*(1. + 2.*rr/rho2));

    ggcov[2][2] = rho2;

    ggcov[3][0] = ggcov[0][3] ;
    ggcov[3][1]  = ggcov[1][3] ;
    ggcov[3][3] = s2*(rho2 + a*a*s2*(1. + 2.*rr/rho2));

    //convert to code coordinates
    for(i=0;i<NDIM;i++){
      for(j=0;j<NDIM;j++){
        gcovp[i][j] = 0;
        for(k=0;k<NDIM;k++) {
          for(l=0;l<NDIM;l++){
            gcovp[i][j] += ggcov[k][l]*ddxdxp[k][i]*ddxdxp[l][j];
          }
        }
      }
    }
  }
  else {
    gcovp[0][0] = -1. ;
    gcovp[1][1] = 1. ;
    gcovp[2][2] = 1. ;
    gcovp[3][3] = 1. ;
  }
}

void gcon_func(double gcov[][NDIM], double gcon[][NDIM])
{
  invert_matrix( gcov, gcon );
}

double gdet_func(double gcov[][NDIM])
{
  int i,j,k;
  int permute[NDIM];
  double gcovtmp[NDIM][NDIM];
  double detg;

  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  gcovtmp[0][i] = gcov[0][i]; }
  if( LU_decompose( gcovtmp,  permute ) != 0  ) {
    fprintf(stderr, "gdet_func(): singular matrix encountered! \n");
  }
  detg = 1.;
  DLOOPA detg *= gcovtmp[j][j];
  return( sqrt(fabs(detg)) );

}

void get_connection(double *X, struct of_geom *geom, double conn[][NDIM][NDIM])
{
	int i,j,k,l ;
	double tmp[NDIM][NDIM][NDIM] ;
	double Xh[NDIM],Xl[NDIM] ;
	double gh[NDIM][NDIM] ;
	double gl[NDIM][NDIM] ;

	for(k=0;k<NDIM;k++) {
		for(l=0;l<NDIM;l++) Xh[l] = X[l] ;
		for(l=0;l<NDIM;l++) Xl[l] = X[l] ;
		Xh[k] += DELTA ;
		Xl[k] -= DELTA ;
		gcov_func(Xh,gh) ;
		gcov_func(Xl,gl) ;

		for(i=0;i<NDIM;i++)
		for(j=0;j<NDIM;j++)
			conn[i][j][k] = (gh[i][j] - gl[i][j])/(Xh[k] - Xl[k]) ;
	}

	// now rearrange to find \Gamma_{ijk}
	for(i=0;i<NDIM;i++)
	for(j=0;j<NDIM;j++)
	for(k=0;k<NDIM;k++)
		tmp[i][j][k] = 0.5*(conn[j][i][k] + conn[k][i][j] - conn[k][j][i]) ;

	// finally, raise index
	for(i=0;i<NDIM;i++)
	for(j=0;j<NDIM;j++)
	for(k=0;k<NDIM;k++)  {
		conn[i][j][k] = 0. ;
		for(l=0;l<NDIM;l++) conn[i][j][k] += geom->gcon[i][l]*tmp[l][j][k] ;
	}
}

void coord(int i, int j, int k, double *X)
{
    // Return zone-centered values for coordinates

	X[0] = startx[0];
	X[1] = startx[1] + (i + 0.5) * dx[1];
	X[2] = startx[2] + (j + 0.5) * dx[2];
	X[3] = startx[3] + (k + 0.5) * dx[3];
}

void get_geometry(int ii, int jj, int kk, struct of_geom *geom)
{
	int i;

	for (i = 0; i <= NDIM*NDIM - 1; i++) {
        geom->gcon[0][i] = grid_gcon[ii][jj][kk][0][i];
        geom->gcov[0][i] = grid_gcov[ii][jj][kk][0][i];
	}
	geom->g = grid_gdet[ii][jj][kk];
}

void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM])
{
	double pphi = fmod(X[3], stopx[3]);
  	if (pphi < 0.0) pphi = stopx[3]+pphi;

	int abc;

	*i = (int) ((X[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
	*j = (int) ((X[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;
	//*k = (int) ((X[3] - startx[3]) / dx[3] - 0.5 + 1000) - 1000;
	*k = (int) ((pphi - startx[3]) / dx[3] - 0.5 + 1000) - 1000;

	abc = *k;

	if (*i < 0) {
		*i = 0;
		del[1] = 0.;
	} else if (*i > N1 - 2) {
		*i = N1 - 2;
		del[1] = 1.;
	} else {
		del[1] = (X[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
	}

	if (*j < 0) {
		*j = 0;
		del[2] = 0.;
	} else if (*j > N2 - 2) {
		*j = N2 - 2;
		del[2] = 1.;
	} else {
		del[2] = (X[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2];
	}

	if (*k < 0) {
		*k = 0;
		del[3] = 0.;
	} else if (*k > N3 - 2) {
		*k = N3 - 2;
		del[3] = 1.;
	} else {
		del[3] = (pphi - ((*k + 0.5) * dx[3] + startx[3])) / dx[3];
	}

	//fprintf(stderr, "Xtoijk: %g %g %g\n", del[1], del[2], del[3]);

	if (del[1] > 1.) del[1] = 1.;
  	if (del[1] < 0.) del[1] = 0.;
  	if (del[2] > 1.) del[2] = 1.;
  	if (del[2] < 0.) del[2] = 0.;
  	if (del[3] > 1.) del[3] = 1.;
  	if (del[3] < 0.) {
    	int oldk = *k;
    	*k = N3-1;
    	del[3] += 1.;
    	if (del[3] < 0) {
      		fprintf(stderr, "Xtoijk: unable to resolve X[3] coordinate to zone %d %d %g %g\n", oldk, *k, del[3], X[3]);
      		exit(-7);
    	}
  	}
	//printf("%g %d %d %g\n", phi, abc, *k, del[3]);


	return;
}


double dot3(double *v1, double *v2)
{
	double vv1[4], vv2[4];

	vv1[0] = 0.;
	vv1[1] = v1[0];
	vv1[2] = v1[1];
	vv1[3] = v1[2];

	vv2[0] = 0.;
	vv2[1] = v2[0];
	vv2[2] = v2[1];
	vv2[3] = v2[2];

	return vv1[0]*vv2[0] + vv1[1]*vv2[1] + vv1[2]*vv2[2] + vv1[3]*vv2[3];
}


void bl_coord(double *X, double *r, double *th, double *phi)
{
  double V[NDIM];
  bl_coord_vec(X,V);
  *r = V[1];
  *th = V[2];
  *phi = V[3];
  return ;
}


/****************************************/

void dxdxp_func(double *X, double dxdxp[][NDIM])
{
  int i,j,k,l ;
  double Xh[NDIM],Xl[NDIM] ;
  double Vh[NDIM],Vl[NDIM] ;

  if(BL){
    for(k=0;k<NDIM;k++) {
      for(l=0;l<NDIM;l++) Xh[l] = X[l] ;
      for(l=0;l<NDIM;l++) Xl[l] = X[l] ;
      Xh[k] += DELTA ;
      Xl[k] -= DELTA ;

      bl_coord_vec(Xh,Vh) ;
      bl_coord_vec(Xl,Vl) ;

      for(j=0;j<NDIM;j++){
        dxdxp[j][k] = (Vh[j]-Vl[j])/(Xh[k] - Xl[k]) ;
        //fprintf(stderr, "%lf\n", dxdxp[j][k]);
      }
    }
  }
  else{
    for(i=0;i<NDIM;i++) {
      for(j=0;j<NDIM;j++) {
        dxdxp[i][j] = 0.;
      }
    }
  }
}

/******************************************************************************/

int invert_matrix( double Am[][NDIM], double Aminv[][NDIM] )
{
  // Uses LU decomposition and back substitution to invert a matrix
  // A[][] and assigns the inverse to Ainv[][]. This routine does not
  // destroy the original matrix A[][].
  // Returns (1) if a singular matrix is found,  (0) otherwise.

  int i,j;
  int n = NDIM;
  int permute[NDIM];
  double dxm[NDIM], Amtmp[NDIM][NDIM];

  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  Amtmp[0][i] = Am[0][i]; }

  // Get the LU matrix:
  if( LU_decompose( Amtmp,  permute ) != 0  ) {
    fprintf(stderr, "invert_matrix(): singular matrix encountered! \n");
    return(1);
  }

  for( i = 0; i < n; i++ ) {
    for( j = 0 ; j < n ; j++ ) { dxm[j] = 0. ; }
    dxm[i] = 1.;

    /* Solve the linear system for the i^th column of the inverse matrix: :  */
    LU_substitution( Amtmp,  dxm, permute );

    for( j = 0 ; j < n ; j++ ) {  Aminv[j][i] = dxm[j]; }

  }

  return(0);
}


int LU_decompose( double A[][NDIM], int permute[] )
{
  // Performs a LU decomposition of the matrix A using Crout's method
  // with partial implicit pivoting.  The exact LU decomposition of the
  // matrix can be reconstructed from the resultant row-permuted form via
  // the integer array permute[]
  // The algorithm closely follows ludcmp.c of "Numerical Recipes
  // in C" by Press et al. 1992.
  // This will be used to solve the linear system  A.x = B
  // Returns (1) if a singular matrix is found,  (0) otherwise.


  const  double absmin = 1.e-30; /* Value used instead of 0 for singular matrices */

  static double row_norm[NDIM];
  double  absmax, maxtemp, mintemp;

  int i, j, k, max_row;
  int n = NDIM;


  max_row = 0;

  // Find the maximum elements per row so that we can pretend later
  // we have unit-normalized each equation:

  for( i = 0; i < n; i++ ) {
    absmax = 0.;

    for( j = 0; j < n ; j++ ) {

      maxtemp = fabs( A[i][j] );

      if( maxtemp > absmax ) {
	absmax = maxtemp;
      }
    }

    // Make sure that there is at least one non-zero element in this row:
    if( absmax == 0. ) {
     fprintf(stderr, "LU_decompose(): row-wise singular matrix!\n");
      return(1);
    }

    row_norm[i] = 1. / absmax ;   // Set the row's normalization factor.
  }

  // The following the calculates the matrix composed of the sum
  // of the lower (L) tridagonal matrix and the upper (U) tridagonal
  // matrix that, when multiplied, form the original maxtrix.
  // This is what we call the LU decomposition of the maxtrix.
  // It does this by a recursive procedure, starting from the
  // upper-left, proceding down the column, and then to the next
  // column to the right.  The decomposition can be done in place
  // since element {i,j} require only those elements with {<=i,<=j}
  // which have already been computed.
  // See pg. 43-46 of "Num. Rec." for a more thorough description.


  // For each of the columns, starting from the left ...
  for( j = 0; j < n; j++ ) {

    // For each of the rows starting from the top....

    // Calculate the Upper part of the matrix:  i < j :
    for( i = 0; i < j; i++ ) {
      for( k = 0; k < i; k++ ) {
	A[i][j] -= A[i][k] * A[k][j];
      }
    }

    absmax = 0.0;

    // Calculate the Lower part of the matrix:  i <= j :

    for( i = j; i < n; i++ ) {

      for (k = 0; k < j; k++) {
	A[i][j] -= A[i][k] * A[k][j];
      }

      // Find the maximum element in the column given the implicit
	  // unit-normalization (represented by row_norm[i]) of each row:
      maxtemp = fabs(A[i][j]) * row_norm[i] ;

      if( maxtemp >= absmax ) {
	absmax = maxtemp;
	max_row = i;
      }

    }

    // Swap the row with the largest element (of column j) with row_j.  absmax
    // This is the partial pivoting procedure that ensures we don't divide
    // by 0 (or a small number) when we solve the linear system.
    // Also, since the procedure starts from left-right/top-bottom,
    // the pivot values are chosen from a pool involving all the elements
    // of column_j  in rows beneath row_j.  This ensures that
    // a row  is not permuted twice, which would mess things up.
    if( max_row != j ) {

      // Don't swap if it will send a 0 to the last diagonal position.
	  // Note that the last column cannot pivot with any other row,
	  // so this is the last chance to ensure that the last two
	  // columns have non-zero diagonal elements.

      if( (j == (n-2)) && (A[j][j+1] == 0.) ) {
	max_row = j;
      }
      else {
	for( k = 0; k < n; k++ ) {

	  maxtemp       = A[   j   ][k] ;
	  A[   j   ][k] = A[max_row][k] ;
	  A[max_row][k] = maxtemp;

	}

	/* Don't forget to swap the normalization factors, too...
	   but we don't need the jth element any longer since we
	   only look at rows beneath j from here on out.
	*/
	row_norm[max_row] = row_norm[j] ;
      }
    }

    /* Set the permutation record s.t. the j^th element equals the
       index of the row swapped with the j^th row.  Note that since
       this is being done in successive columns, the permutation
       vector records the successive permutations and therefore
       index of permute[] also indexes the chronology of the
       permutations.  E.g. permute[2] = {2,1} is an identity
       permutation, which cannot happen here though.
    */

    permute[j] = max_row;

    if( A[j][j] == 0. ) {
      A[j][j] = absmin;
    }


  /* Normalize the columns of the Lower tridiagonal part by their respective
     diagonal element.  This is not done in the Upper part because the
     Lower part's diagonal elements were set to 1, which can be done w/o
     any loss of generality.
  */
    if( j != (n-1) ) {
      maxtemp = 1. / A[j][j]  ;

      for( i = (j+1) ; i < n; i++ ) {
	A[i][j] *= maxtemp;
      }
    }

  }

  return(0);

  /* End of LU_decompose() */

}


/************************************************************************

   LU_substitution():

       Performs the forward (w/ the Lower) and backward (w/ the Upper)
       substitutions using the LU-decomposed matrix A[][] of the original
       matrix A' of the linear equation:  A'.x = B.  Upon entry, A[][]
       is the LU matrix, B[] is the source vector, and permute[] is the
       array containing order of permutations taken to the rows of the LU
       matrix.  See LU_decompose() for further details.

       Upon exit, B[] contains the solution x[], A[][] is left unchanged.

************************************************************************/

void LU_substitution( double A[][NDIM], double B[], int permute[] )
{
  int i, j ;
  int n = NDIM;
  double tmpvar,tmpvar2;


  /* Perform the forward substitution using the LU matrix.
   */
  for(i = 0; i < n; i++) {

    /* Before doing the substitution, we must first permute the
       B vector to match the permutation of the LU matrix.
       Since only the rows above the currrent one matter for
       this row, we can permute one at a time.
    */
    tmpvar        = B[permute[i]];
    B[permute[i]] = B[    i     ];
    for( j = (i-1); j >= 0 ; j-- ) {
      tmpvar -=  A[i][j] * B[j];
    }
    B[i] = tmpvar;
  }


  /* Perform the backward substitution using the LU matrix.
   */
  for( i = (n-1); i >= 0; i-- ) {
    for( j = (i+1); j < n ; j++ ) {
      B[i] -=  A[i][j] * B[j];
    }
    B[i] /= A[i][i] ;
  }

  /* End of LU_substitution() */

}


/****************************************/

void bl_coord_vec(double *X, double *V)
{
    void vofx_gammiecoords(double *X, double *V);
    void vofx_sjetcoords( double *X, double *V );

    if (BL == 1) {
        if(!DOCYLINDRIFYCOORDS) {
            //use original coordinates
            vofx_gammiecoords( X, V );
        }
        else {
            //apply cylindrification to original coordinates
            //[this internally calls vofx_gammiecoords()]
            vofx_cylindrified( X, vofx_gammiecoords, V );
        }
    }
    /*
    if (BL == 1) {
#if(!DOCYLINDRIFYCOORDS)
	    //use original coordinates
        vofx_gammiecoords( X, V );
        printf("GAMMIE\n");
#else
        //apply cylindrification to original coordinates
        //[this internally calls vofx_gammiecoords()]
        vofx_cylindrified( X, vofx_gammiecoords, V );
        printf("CYL\n");
#endif
}*/
	// NOT THE CASE IN THIS INITIAL VERSION! -- GS
//	else if (BL == 2) {
//#if(!DOCYLINDRIFYCOORDS)
        //use original coordinates
//        vofx_sjetcoords( X, V );
//#else
        //apply cylindrification to original coordinates
        //[this internally calls vofx_sjetcoords()]
//        vofx_cylindrified( X, vofx_sjetcoords, V );
//#endif
    //}
    else {
        V[0] = X[0];
        V[1] = X[1];
        V[2] = X[2];
        V[3] = X[3];
    }
	// avoid singularity at polar axis
#if(COORDSINGFIX && BL)
    if(fabs(V[2])<SINGSMALL) {
  		if((V[2])>=0) V[2] =  SINGSMALL;
  		if((V[2])<0)  V[2] = -SINGSMALL;
	}
	if(fabs(M_PI - V[2]) < SINGSMALL) {
	  	if((V[2])>=M_PI) V[2] = M_PI+SINGSMALL;
	  	if((V[2])<M_PI)  V[2] = M_PI-SINGSMALL;
	}
#endif

	return;
}

void vofx_gammiecoords(double *X, double *V)
{
	// input: code coordinates (MKS) X
	// output: KS/BL coordinates V

	double theexp = X[1];

	if (X[1] > x1br) {
		// hyperexponential for X[1] > x1br
    	theexp += cpow2 * pow(X[1] - x1br, npow2);
  	}

	// calculate KS/BL coordinates based on MKS input X
	V[0] = X[0];
	V[1] = exp(theexp) + R0;
	V[2] = M_PI_2*(1.0 + X[2]) + ((1. - hslope)/2.)*sin(M_PI*(1.0 + X[2]));
	V[3] = X[3];
}

void vofx_sjetcoords( double *X, double *V )
{
    //for SJETCOORDS
    double theexp;
    double Ftrgen( double x, double xa, double xb, double ya, double yb );
    double limlin( double x, double x0, double dx, double y0 );
    double minlin( double x, double x0, double dx, double y0 );
    double mins( double f1, double f2, double df );
    double maxs( double f1, double f2, double df );
    double thetaofx2(double x2, double ror0nu);
    double  fac, faker, ror0nu;
    double fakerdisk, fakerjet;
    double rbeforedisk, rinsidedisk, rinsidediskmax, rafterdisk;
    double ror0nudisk, ror0nujet, thetadisk, thetajet;

    V[0] = X[0];

    theexp = X[1];

    if( X[1] > x1br ) {
        theexp += cpow2 * pow(X[1]-x1br,npow2);
    }
    V[1] = R0+exp(theexp);

#if(0) //JON's method
    myhslope=2.0-Qjet*pow(V[1]/r1jet,-njet*(0.5+1.0/M_PI*atan(V[1]/r0grid-rsjet/r0grid)));

    if (X[2] < 0.5) {
        V[2] = M_PI * X[2] + ((1. - myhslope) / 2.) * mysin(2. * M_PI * X[2]);
    }
    else {
        V[2] = M_PI - (M_PI * (1.0-X[2])) + ((1. - myhslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
    }
#elif(1) //SASHA's
    fac = Ftrgen( fabs(X[2]), fracdisk, 1-fracjet, 0, 1 );

    rbeforedisk = mins( V[1], r0disk, 0.5*r0disk );
    rinsidedisk = 1.;
    rinsidediskmax = 1.;
    rafterdisk = maxs( 1, 1 + (V[1]-rdiskend)*r0jet/(rjetend*r0disk*rinsidediskmax), 0.5*rdiskend*r0jet/(rjetend*r0disk*rinsidediskmax) );

    fakerdisk = rbeforedisk * rinsidedisk * rafterdisk;

    fakerjet = mins( V[1], r0jet, 0.5*r0jet ) * maxs( 1, V[1]/rjetend, 0.5 );

    ror0nudisk = pow( (fakerdisk - rsjet*Rin)/r0grid, jetnu/2 );
    ror0nujet = pow( (fakerjet - rsjet*Rin)/r0grid, jetnu/2 );
    thetadisk = thetaofx2( X[2], ror0nudisk );
    thetajet = thetaofx2( X[2], ror0nujet );
    V[2] = fac*thetajet + (1 - fac)*thetadisk;
#else
    V[2] = M_PI_2 * (1.0 + X[2]);
#endif

    // default is uniform \phi grid
    V[3]=X[3];
}

double thetaofx2(double x2, double ror0nu)
{
    double theta;
    if( x2 < -0.5 ) {
        theta = 0       + atan( tan((x2+1)*M_PI_2)/ror0nu );
    }
    else if( x2 >  0.5 ) {
        theta = M_PI    + atan( tan((x2-1)*M_PI_2)/ror0nu );
    }
    else {
        theta = M_PI_2 + atan( tan(x2*M_PI_2)*ror0nu );
    }
    return(theta);
}

/*** CYLINDRIFICATION ROUTINES ***/

double Ftr( double x )
{
    double res;

    if ( x <= 0 ) {
        res = 0;
    }
    else if( x >= 1 ) {
        res = 1;
    }
    else {
        res = (64 + cos(5*M_PI*x) + 70*sin((M_PI*(-1 + 2*x))/2.) + 5*sin((3*M_PI*(-1 + 2*x))/2.))/128.;
    }

    return( res );
}

double Ftrgenlin( double x, double xa, double xb, double ya, double yb )
{
    double Ftr( double x );
    double res;

    res = (x*ya)/xa + (-((x*ya)/xa) + ((x - xb)*(1 - yb))/(1 - xb) + yb)*Ftr((x - xa)/(-xa + xb));

    return( res );
}

//goes from ya to yb as x goes from xa to xb
double Ftrgen( double x, double xa, double xb, double ya, double yb )
{
    double Ftr( double x );
    double res;

    res = ya + (yb-ya)*Ftr( (x-xa)/(xb-xa) );

    return( res );
}

double Fangle(double x) {
	// Equation B3 of Ressler+2017
	// input: angle in MKS coordinates
	// output: "smoothified" angle

	double res;

  	if(x <= -1) {
    	res = 0;
  	}
  	else if (x >= 1) {
    	res = x;
  	}
  	else {
    	res = (1 + x + (-140*sin((M_PI*(1 + x))/2.) +
            (10*sin((3*M_PI*(1 + x))/2.))/3. +
            (2*sin((5*M_PI*(1 + x))/2.))/5.)/(64.*M_PI))/2.;
  	}
    return(res);
}



double limlin(double x, double x0, double dx, double y0) {
	// RHS of Equation B4 of Ressler+2017.

  	double Fangle(double x);
  	return(y0 - dx * Fangle(-(x-x0)/dx));
}

double mins(double f1, double f2, double df) {
	// Equation B4 in Ressler+2017

    double limlin(double x, double x0, double dx, double y0);
    return(limlin(f1, f2, df, f2));
}

double maxs(double f1, double f2, double df) {
	// Equation B5 in Ressler+2017

    double mins(double f1, double f2, double df);
    return(-mins(-f1, -f2, df));
}

//=mins if dir < 0
//=maxs if dir >= 0
double minmaxs(double f1, double f2, double df, double dir) {
	// Returns either min_s (B4 in Ressler+2017) or max_s (B5 in Ressler+2017)

    double mins(double f1, double f2, double df);
    double maxs(double f1, double f2, double df);
    if (dir>=0) {
    	return(maxs(f1, f2, df));
    }
  	return(mins(f1, f2, df));
}

static double sinth0(double *X0, double *X, void (*vofx)(double*, double*));
static double sinth1in(double *X0, double *X, void (*vofx)(double*, double*));
static double th2in(double *X0, double *X, void (*vofx)(double*, double*));
static void to1stquadrant(double *Xin, double *Xout, int *ismirrored);
static double func1(double *X0, double *X,  void (*vofx)(double*, double*));
static double func2(double *X0, double *X,  void (*vofx)(double*, double*));

//Converts copies Xin to Xout and converts
//but sets Xout[2] to lie in the 1st quadrant, i.e. Xout[2] \in [-1,0])
//if the point had to be mirrored
void to1stquadrant(double *Xin, double *Xout, int *ismirrored) {
	// input: coordinates in MKS (Xin)
	// output: coordinates in MKS (Xout) lying in the 1st quadrant

  	double ntimes;
	int j;

  	DLOOPA Xout[j] = Xin[j];

  	//bring the angle variables to -2..2 (for X) and -2pi..2pi (for V)
  	ntimes = floor((Xin[2]+2.0)/4.0);
  	//this forces -2 < Xout[2] < 2
  	Xout[2] -= 4 * ntimes;

  	*ismirrored = 0;

  	if( Xout[2] > 0. ) {
    	Xout[2] = -Xout[2];
    	*ismirrored = 1-*ismirrored;
  	}

  	//now force -1 < Xout[2] < 0
  	if( Xout[2] < -1. ) {
    	Xout[2] = -2. - Xout[2];
    	*ismirrored = 1-*ismirrored;
  	}
}

double sinth0(double *X0, double *X, void (*vofx)(double*, double*)) {
	// Used in th0 = asin(sinth0(X0, X, vofx));
	// input: X0 (MKS cyl), X (MKS mirrored)
	// output: sine of smoothened angle

	double V0[NDIM];
	double Xc0[NDIM], Vc0[NDIM];
    int j;

    //X1 = {0, X[1], X0[1], 0}
	// Copy input X (mirrored MKS) into Xc0
    DLOOPA Xc0[j] = X[j];
	// Set Xc0 to be cylindrified MKS
    Xc0[2] = X0[2];

	// Finds associated Vc0 (KS/BL) to Xc0 (MKS mirrored, but Xc0[2] cylindrified)
    vofx(Xc0, Vc0);
	// Finds associated V0 (KS/BL) to X0 (MKS cylindrified)
    vofx(X0, V0);

	//    (r cyl)* sin(theta cyl)/r mirrored
    return(V0[1] * sin(V0[2]) / Vc0[1]);
}

double sinth1in(double *X0, double *X, void (*vofx)(double*, double*)) {
	// input:
	// output:

	double V0[NDIM];
	double V[NDIM];
    double X0c[NDIM], V0c[NDIM];
    int j;

    //X1 = {0, X[1], X0[1], 0}
    DLOOPA X0c[j] = X0[j];
    X0c[2] = X[2];

    vofx(X, V);
    vofx(X0c, V0c);
    vofx(X0, V0);

    return(V0[1] * sin(V0c[2]) / V[1]);
}


double th2in(double *X0, double *X, void (*vofx)(double*, double*)) {
    double V0[NDIM];
    double V[NDIM];
    double Xc0[NDIM], Vc0[NDIM];
    double Xcmid[NDIM], Vcmid[NDIM];
    int j;
    double res, th0;

	// Copy input X (mirrored MKS) into Xc0
    DLOOPA Xc0[j] = X[j];
	// Define Xc0[2] as X0[2] (cylindrified, global_x20)
    Xc0[2] = X0[2];
	// Find the associated Vc0 (KS/BL) to Xc0
    vofx(Xc0, Vc0); // vofx_gammiecoords in our case
	// Xc0 = (X[0], X[1], global_x20, X[3])
	// Vc0 = KS/BL of Xc0
	// So Vc0[2] = M_PI_2*(1.0+global_x20) + ((1. - hslope)/2.)*sin(M_PI*(1.0+global_x20));

	// Copy input X (mirrored MKS) into Xcmid
    DLOOPA Xcmid[j] = X[j];
	// Set Xcmid[2] to 0
    Xcmid[2] = 0;
	// Find the associated Vcmid (KS/BL) to Xcmid
    vofx(Xcmid, Vcmid); // vofx_gammiecoords in our case

	// Find associated V0 and V to X0 (cylindrified MKS) and X (mirrored MKS)
    vofx(X0, V0); // vofx_gammiecoords in our case
    vofx(X, V); // vofx_gammiecoords in our case

	// Find smoothened angle
    th0 = asin(sinth0(X0, X, vofx));
    res = (V[2] - Vc0[2])/(Vcmid[2] - Vc0[2]) * (Vcmid[2]-th0) + th0;

    return(res);
}

//Adjusts V[2]=theta so that a few innermost cells around the pole
//become cylindrical
//ASSUMES: poles are at
//            X[2] = -1 and +1, which correspond to
//            V[2] = 0 and pi
void vofx_cylindrified(double *Xin, void (*vofx)(double*, double*), double *Vout) {
	// this is grtrans calcth_cylindrified
	// grtrans  ------------ here
	//
	// r0       ------------ V0[1]
	// rin      ------------ V[1]
	// x1tr     ------------ Xtr[1]
	// x20      ------------ X0[2]
	// rtr      ------------ Vtr[1]
	// theta    ------------ Vout[2]
	// thorig   ------------ Vin[2]
	// thmirror ------------ V[2]
	// x2mirror ------------ X[2]
	// calcthmksbl3 -------- vofx
    double npiovertwos;
	double Vin[NDIM];            // KS/BL from input pure MKS Xin
	double X[NDIM], V[NDIM];     // mirrored pure MKS, BL/KS
    double X0[NDIM], V0[NDIM];   // cylindrified MKS, BL/KS
    double Xtr[NDIM], Vtr[NDIM]; // transition MKS, BL/KS
    double f1, f2, dftr;
    double sinth, th;
    int j, ismirrored;


	// In our case, vofx is vofx_gammiecoords, so here we take the input MKS
    // coordinates X and find KS/BL Vin as in vofx_gammiecoords
    vofx(Xin, Vin); // vofx_gammiecoords in our case

    // BRING INPUT TO 1ST QUADRANT:  X[2] \in [-1 and 0]
    to1stquadrant(Xin, X, &ismirrored);
	// After bringing to 1st quadrant (Xin to X -- MKS), we find the associated
	// KS/BL V
    vofx(X, V); // vofx_gammiecoords in our case

    // initialize X0: cylindrify region
    // X[1] < X0[1] && X[2] < X0[2] (value of X0[3] not used)
    X0[0] = Xin[0];
    X0[1] = global_x10;
    X0[2] = global_x20;
    X0[3] = 0;
	// After initializing cylindrification (cylindrified coords being X0 -- MKS)
	// we find the associated KS/BL V0
    vofx(X0, V0); // vofx_gammiecoords in our case

    // {0, roughly midpoint between grid origin and x10, -1, 0}
	// Here, we simply copy the 1st quadrant coordinates X into Xtr
    DLOOPA Xtr[j] = X[j];
	// makes sure that Xtr[1] is always bound to be between startx[1] and X0[1]
    Xtr[1] = log(0.5*( exp(X0[1])+exp(startx[1])));
    // After setting the Xtr coordinates, find associated Vtr
    vofx(Xtr, Vtr); // vofx_gammiecoords in our case

	// calculates V from X and returns sin(V[2]) -- X0 not used?!
    f1 = func1(X0, X, vofx); // vofx_gammiecoords in our case
    //
    f2 = func2(X0, X, vofx); // vofx_gammiecoords in our case
    dftr = func2(X0, Xtr, vofx) - func1(X0, Xtr, vofx); // vofx_gammiecoords in our case

    // Compute new theta
    sinth = maxs(V[1]*f1, V[1]*f2, Vtr[1]*fabs(dftr)+SMALL) / V[1];
    th = asin(sinth);

    //initialize Vout with the original values
    DLOOPA Vout[j] = Vin[j];

    //apply change in theta in the original quadrant
    if (0 == ismirrored) {
	    Vout[2] = Vin[2] + (th - V[2]);
    }
    else {
    //if mirrrored, flip the sign
    	Vout[2] = Vin[2] - (th - V[2]);
	}
}

double func1(double *X0, double *X,  void (*vofx)(double*, double*)) {
	// input:
	// output: returns the sine of MKS-calculated theta

    double V[NDIM];

    vofx(X, V); // vofx_gammiecoords in our case

	return(sin(V[2]));
}

double func2(double *X0, double *X,  void (*vofx)(double*, double*)) {
	// input:
	// output:

    double V[NDIM];
    double Xca[NDIM]; // axis
    double func2;
    int j;
    double sth1in, sth2in, sth1inaxis, sth2inaxis;

    //{0, X[1], -1, 0}
    DLOOPA Xca[j] = X[j];
    Xca[2] = -1;

    vofx(X, V); // vofx_gammiecoords in our case
    // X and V are mirrored pure MKS
    // X0 are cylindrified MKS
    sth1in = sinth1in(X0, X, vofx); // sth1in depends on V[1]
    sth2in = sin(th2in(X0, X, vofx)); // th2in depends on V[2]

    sth1inaxis = sinth1in(X0, Xca, vofx);
    sth2inaxis = sin(th2in(X0, Xca, vofx));

    func2 = minmaxs(sth1in, sth2in, fabs(sth2inaxis-sth1inaxis)+SMALL, X[1] - X0[1]);

    return(func2);
}

/*** END OF CYLINDRIFICATION ROUTINES ***/


#define ITMAX 100
#define EPS 3.0e-8
double zbrent(double (*func)(double, double), double param1, double lower, double upper, double tol)
{
    // Brent method for finding the root of a function, which is known to be
    // between lower and upper. Here, param1 is an additional parameter
    // associated with func.

    // Taken essentially from Numerical Recipes, with few modifications.

	int iter;
	double aa = lower, bb = upper, cc, dd, ee, min1, min2;
	double fa = (*func)(aa, param1),fb = (*func)(bb, param1), fc, pp, qq, rr, ss, tol1, xm;

	//if (fb*fa > 0.0) printf("Root must be bracketed in ZBRENT\n");
	fc = fb;
	for (iter = 1; iter <= ITMAX; iter++) {
		if (fb*fc > 0.0) {
			cc = aa;
			fc = fa;
			ee = dd = bb - aa;
		}
		if (fabs(fc) < fabs(fb)) {
			aa = bb;
			bb = cc;
			cc = aa;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		tol1 = 2.0*EPS*fabs(bb) + 0.5*tol;
		xm = 0.5*(cc - bb);
		if (fabs(xm) <= tol1 || fb == 0.0) return bb;
		if (fabs(ee) >= tol1 && fabs(fa) > fabs(fb)) {
			ss = fb/fa;
			if (aa == cc) {
				pp = 2.0*xm*ss;
				qq = 1.0 - ss;
			} else {
				qq = fa/fc;
				rr = fb/fc;
				pp = ss*(2.0*xm*qq*(qq - rr) - (bb - aa)*(rr - 1.0));
				qq = (qq - 1.0)*(rr - 1.0)*(ss - 1.0);
			}
			if (pp > 0.0)  qq = -qq;
			pp   = fabs(pp);
			min1 = 3.0*xm*qq - fabs(tol1*qq);
			min2 = fabs(ee*qq);
			if (2.0*pp < (min1 < min2 ? min1 : min2)) {
				ee = dd;
				dd = pp/qq;
			} else {
				dd = xm;
				ee = dd;
			}
		} else {
			dd = xm;
			ee = dd;
		}
		aa = bb;
		fa = fb;
		if (fabs(dd) > tol1)
			bb += dd;
		else
			bb += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
		fb = (*func)(bb, param1);
	}
	printf("Maximum number of iterations exceeded in ZBRENT");
}


double find_x1_cyl(double x1, double radius)
{
    // to be used with zbrent
    double theexp = x1;
    if (x1 > x1br) theexp += cpow2 * pow(x1 - x1br, npow2);
    double diff = radius - (exp(theexp) + R0);
    return(diff);
}


double find_x2_cyl(double x2, double theta)
{
    // to be used with zbrent
    double v2 = calcth_cylindrified(x2);
    double diff = theta - v2;
    return(diff);
}


double calcrmks(double x1) {
    double theexp = x1;
    if (x1 > x1br) {
        theexp += cpow2 * pow(x1 - x1br, npow2);
    }
    return(exp(theexp) + R0);
}


double to1stquadrant_single(double x2in)
{
    double ntimes = floor((x2in + 2.)/4.);
    double x2mirror -= 4*ntimes;
    int *ismirrored = 0;

  	if(x2mirror > 0.) {
        x2mirror = -x2mirror;
    	*ismirrored = 1-*ismirrored;
    }

    if (x2mirror < -1.) {
        x2mirror = -2. - x2mirror;
    	*ismirrored = 1-*ismirrored;
    }
    return(x2mirror);
}


double func2_single(double r0, double rr, double x20, double x2)
{
    double mone = -1.;
    double sth1in = sinth1in_single(r0,rr,x20,x2);
    double sth2in = sin(th2in_single(r0,rr,x20,x2));
    double sth1inaxis = sinth1in_single(r0,rr,x20,mone);
    double sth2inaxis = sin(th2in_single(r0,rr,x20,mone));

    return(minmaxs(sth1in, sth2in, abs(sth2inaxis-sth1inaxis)+SMALL, r-r0));
}


double sinth1in_single(double r0, double rr, double x20, double x2)
{
    double thc = M_PI_2*(1.0 + x2) + ((1. - hslope)/2.)*sin(M_PI*(1.0 + x2));
    return (r0 * sin(thc) / rr);
}


double th2in_single(double r0, double rr, double x20, double x2)
{
    double thetac = M_PI_2*(1.0 + x20) + ((1. - hslope)/2.)*sin(M_PI*(1.0 + x20));
    double thetamid = M_PI_2*(1.0 + 0.) + ((1. - hslope)/2.)*sin(M_PI*(1.0 + 0.));
    double ttheta = M_PI_2*(1.0 + x2) + ((1. - hslope)/2.)*sin(M_PI*(1.0 + x2));
    double th0 = asin(sinth0_single(x20, r0, rr));

    return((ttheta - thetac)/(thetamid - thetac) * (thetamid-th0) + th0);
}


double sinth0_single(double x20, double r0, double rr)
{
    double th0 = M_PI_2*(1.0 + x20) + ((1. - hslope)/2.)*sin(M_PI*(1.0 + x20));
    return (r0*sin(th0)/rr);
}


double calcth_cylindrified(double x2in)
{
    double thorig = M_PI_2*(1.0 + x2in) + ((1. - hslope)/2.)*sin(M_PI*(1.0 + x2in));
    double x2mirror = to1stquadrant_single(x2in);
    double thmirror = M_PI_2*(1.0 + x2mirror) + ((1. - hslope)/2.)*sin(M_PI*(1.0 + x2mirror));

    double r0 = calcrmks(global_x10);
    double x20 = global_x20;
    double x1tr = log(0.5*(exp(global_x10)+exp(startx1)));
    double rtr = calcrmks(x1tr);
    double thtr = M_PI_2*(1.0 + x2mirror) + ((1. - hslope)/2.)*sin(M_PI*(1.0 + x2mirror));

    double f1 = sin(M_PI_2*(1.0 + x2mirror) + ((1. - hslope)/2.)*sin(M_PI*(1.0 + x2mirror)));
    double f2 = func2_single(r0, rin, x20, x2mirror);

    double dftr = func2_single(r0, rtr, x20, x2mirror) - f1;
    double sinth = maxs(rin*f1, rin*f2, rtr*abs(dftr)+SMALL)/rin;

    double tth = asin(sinth);
    double ttheta = thorig;

    if (0 == ismirrored) ttheta = thorig + (tth - thmirror);
    else ttheta  = thorig  - (tth  - thmirror);

    return ttheta;
}





/*


double x1_ks_to_x1_mks (double *X_ks)
{
    // finds x1 mks given x1 ks

    double find_x1_mks(double x1, double radius);
    double x1_mks;

    x1_mks = zbrent(find_x1_mks, X_ks[1], 0, 8., 1E-6);

    return x1_mks;
}

double x2_ks_to_x2_mks (double *X_ks)
{
    // finds x2 mks given x2 ks

    double find_x2_mks(double x2, double theta);
    double x2_mks;

    x2_mks = zbrent(find_x2_mks, X_ks[2], -1.0, 1.0, 1E-6);

    return x2_mks;
}

double find_x1_mks(double x1, double radius)
{
    // given ks radius, find mks x1
    // to be used with zbrent

    double diff = radius - (exp(x1) + R0);
    return(diff);
}

double find_x2_mks(double x2, double theta)
{
    // given ks theta, find mks x2
    // to be used with zbrent

    double diff = theta - (M_PI_2*(1.0 + x2) + ((1. - hslope)/2.)*sin(M_PI*(1.0 + x2)));
    return(diff);
}

void x_cyl_to_ks(double *X_cyl, double *X_ks)
{
    // transforms coordinates from cyl to ks

    double X_aux[NDIM];

    bl_coord_vec(X_cyl, X_aux);

    X_ks[0] = X_aux[0];
    X_ks[1] = X_aux[1];
    X_ks[2] = X_aux[2];
    X_ks[3] = X_aux[3];
}

void x_ks_to_mks(double *X_ks, double *X_mks)
{
    // transforms coordinates from ks to mks

    X_mks[0] = X_ks[0];
    X_mks[1] = x1_ks_to_x1_mks(X_ks);
    X_mks[2] = x2_ks_to_x2_mks(X_ks);
    X_mks[3] = X_ks[3];
}



void jac_cyl_to_ks(double *X_cyl, double jac[][NDIM])
{
    // calculates jacobian from cyl to ks
    // ks = d(ks)/d(cyl) * cyl

    double Xp_cyl[NDIM], Xm_cyl[NDIM];
    double Xp_ks[NDIM], Xm_ks[NDIM];
    int i, j;

    // DELTAs, Xp, Xm found after some trial and error
    Xp_cyl[0] = X_cyl[0] + DELTA;
    Xm_cyl[0] = X_cyl[0] - DELTA;
    Xp_cyl[1] = X_cyl[1]*(1. + 0.5*DELTAR);
    Xm_cyl[1] = X_cyl[1]*(1. - 0.5*DELTAR);
    Xp_cyl[2] = X_cyl[2] + DELTA;
    Xm_cyl[2] = X_cyl[2] - DELTA;
    Xp_cyl[3] = X_cyl[3] + DELTA;
    Xm_cyl[3] = X_cyl[3] - DELTA;

    bl_coord_vec(Xp_cyl, Xp_ks);
    bl_coord_vec(Xm_cyl, Xm_ks);

    for (i = 0; i < NDIM; i++) {
        for (j = 0; j < NDIM; j++) {
            jac[i][j] = 0.;
        }
    }

    jac[0][0] = 1.;
    jac[1][1] = (Xp_ks[1] - Xm_ks[1])/(Xp_cyl[1] - Xm_cyl[1]);
    jac[2][1] = (Xp_ks[2] - Xm_ks[2])/(Xp_cyl[1] - Xm_cyl[1]);
    jac[2][2] = (Xp_ks[2] - Xm_ks[2])/(Xp_cyl[2] - Xm_cyl[2]);
    jac[3][3] = 1.;
}

void jac_ks_to_mks(double *X_ks, double jac[][NDIM])
{
    // calculates jacobian from ks to mks
    // mks = d(mks)/d(ks) * ks

    double X_mks[NDIM];
    double Xp_mks[NDIM], Xm_mks[NDIM];
    double Xp_ks[NDIM], Xm_ks[NDIM];
    int i, j;

    for (i = 0; i < NDIM; i++) {
        for (j = 0; j < NDIM; j++) {
            jac[i][j] = 0.;
        }
    }

    X_mks[1] = x1_ks_to_x1_mks(X_ks);
    X_mks[2] = x2_ks_to_x2_mks(X_ks);

    jac[0][0] = 1.;
    jac[1][1] = 1./(exp(X_mks[1]));
	jac[2][2] = 1./(M_PI_2 - 0.5 * M_PI * (hslope - 1.) * cos(M_PI * (X_mks[2] + 1.)));
    jac[3][3] = 1.;
}

void fourvec_old_to_new(double *fourvec_old, double jac[][NDIM], double *fourvec_new)
{
    // given a 4-vector in X coordinates, calculates same 4-vector in Y
    // coordinates using jacobian from X to Y

    int i;

    for  (i = 0; i < NDIM; i++) {
        fourvec_new[i] = 0;
    }

    for (i = 0; i < NDIM; i++) {
        fourvec_new[0] += jac[i][0]*fourvec_old[i];
        fourvec_new[1] += jac[i][1]*fourvec_old[i];
        fourvec_new[2] += jac[i][2]*fourvec_old[i];
        fourvec_new[3] += jac[i][3]*fourvec_old[i];
    }
}

*/
