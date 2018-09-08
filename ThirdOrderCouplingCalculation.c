/*
 * This function finds the mode combinations which lead to non-zero coupling constant
 */

#include <omp.h>
#include "struct_def.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

/* Definitions for routines */
double GalerkinProj(gsl_vector* , gsl_vector*, gsl_vector* , disc_const );
void ModeshapesDisp(gsl_vector* , int , int , disc_const );

void ThirdOrderCouplingCalculation(sys_var* s, mat_const mc, disc_const dc) {
  printf("#    - Third order mode coupling calculations ...\n");

  /**************************/
  // Third order mode
  // coupling!
  /**************************/

  int i, j, k, l, n;
  int cou1, cou2, cou3, cou4;
  int sind, pind, qind, rind;
  int idat;
  FILE *infp, *outfp;

  double fdat, gams_pqr, Hn_pq, Hn_rs, xeta_n4, norm2D;
  int mords, nords, mordp, nordp, mordq, nordq, mordr, nordr;

  gsl_matrix * Fmat 	  = gsl_matrix_alloc (dc.Nfx*dc.Nfy, dc.Nfx*dc.Nfy);
  gsl_vector * Feval 	  = gsl_vector_alloc (dc.Nfx*dc.Nfy);
  gsl_matrix * Fevec 	  = gsl_matrix_alloc (dc.Nfx*dc.Nfy, dc.Nfx*dc.Nfy);
  s->alphamat  	  	  = gsl_matrix_calloc (s->IRcou, 2);
  // for (i = 0; i < s->IRcou; i++){
  //   gsl_matrix_set(s->alphamat, i, 0, 0.0);
  // }

  gsl_vector *phi_s, *phi_p, *phi_q, *phi_r, *psi_n;

  // First compute the eigen vectors and eigen values of the 
  // FD discretized BH operator.
  // Discretization is done using Mathematica
  // Fmat import !!

  printf("#         o Input file required: 1. in.Fmat\n");
  // check for the required input files

  if ((infp=fopen("in.Fmat", "r")) == NULL) 
  {
     printf("#           ERROR: File not found! \n");
     exit(1);
  }
  infp = fopen("in.Fmat", "r");
  for (i = 0; i < dc.Nfx*dc.Nfy; i++) {
    for (j = 0; j < dc.Nfx*dc.Nfy; j++) {
      fscanf(infp, "%lf", &fdat);
      gsl_matrix_set(Fmat, i, j, fdat);
    }
  }
  fclose(infp);
  // outfp = fopen("Fmat.check", "w");
  // for (i = 0; i < dc.Nfx*dc.Nfy; i++) {
  //   for (j = 0; j < dc.Nfx*dc.Nfy; j++) {
  //     fprintf(infp, "%5.5e\t", gsl_matrix_get(Fmat, i, j));
  //   }
  //   fprintf(infp, "\n");
  // }
  // fclose(outfp);

  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (dc.Nfx*dc.Nfy);
  gsl_eigen_symmv (Fmat, Feval, Fevec, w);
  gsl_eigen_symmv_free (w);
  gsl_eigen_symmv_sort (Feval, Fevec, GSL_EIGEN_SORT_ABS_ASC);

  printf("#          Mode groups and couling constants ...\n");
  // First we will set the values of alphamat corresponding to first occurance of ids as 0.0 and all duplicate occurance of ids as -100.0
  for (i = 0; i < s->IRcou; i++){
    idat = (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 4);
    fdat =       gsl_matrix_get(s->alphamat, i, 0);
    if (fdat==0.0){
      for (j = i+1; j < s->IRcou; j++){
        if ((int) gsl_matrix_get(s->IRs_pqrcombmatsorted, j, 4)==idat){
	  gsl_matrix_set(s->alphamat, j, 0, -100.0);
        }
      }
    }
  }

  #pragma omp parallel for private(i,sind,pind,qind,rind,idat,fdat,mords,nords,phi_s,norm2D,mordp,nordp,phi_p,mordq,nordq,phi_q,mordr,nordr,phi_r,gams_pqr,n,psi_n,xeta_n4,Hn_pq,Hn_rs,j)
  for (i = 0; i < s->IRcou; i++){
    // printf("Thread %d is doing iteration %d.\n", omp_get_thread_num( ), i);
    
    phi_s 	  = gsl_vector_alloc(dc.Nfx*dc.Nfy);
    phi_p 	  = gsl_vector_alloc(dc.Nfx*dc.Nfy);
    phi_q 	  = gsl_vector_alloc(dc.Nfx*dc.Nfy);
    phi_r 	  = gsl_vector_alloc(dc.Nfx*dc.Nfy);
    psi_n 	  = gsl_vector_alloc(dc.Nfx*dc.Nfy);
    fdat =       gsl_matrix_get(s->alphamat, i, 0);
    // looping over unique ids
    if (fdat==0.0){ 
      sind = (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 0);
      pind = (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 1);
      qind = (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 2);
      rind = (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 3);
      idat = (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, i, 4);

      mords = gsl_matrix_get(s->modindmat, sind-1, 0); 
      nords = gsl_matrix_get(s->modindmat, sind-1, 1);
      ModeshapesDisp(phi_s, mords, nords, dc);
      norm2D  = gsl_blas_dnrm2(phi_s)*pow(dc.hx*dc.hy, 0.5);
      gsl_vector_scale (phi_s, 1.0/norm2D);
   
      mordp = gsl_matrix_get(s->modindmat, pind-1, 0); 
      nordp = gsl_matrix_get(s->modindmat, pind-1, 1);
      ModeshapesDisp(phi_p, mordp, nordp, dc);
      norm2D  = gsl_blas_dnrm2(phi_p)*pow(dc.hx*dc.hy, 0.5);
      gsl_vector_scale (phi_p, 1.0/norm2D);
  
      mordq = gsl_matrix_get(s->modindmat, qind-1, 0); 
      nordq = gsl_matrix_get(s->modindmat, qind-1, 1);
      ModeshapesDisp(phi_q, mordq, nordq, dc);
      norm2D  = gsl_blas_dnrm2(phi_q)*pow(dc.hx*dc.hy, 0.5);
      gsl_vector_scale (phi_q, 1.0/norm2D);
  
      mordr = gsl_matrix_get(s->modindmat, rind-1, 0); 
      nordr = gsl_matrix_get(s->modindmat, rind-1, 1);
      ModeshapesDisp(phi_r, mordr, nordr, dc);
      norm2D  = gsl_blas_dnrm2(phi_r)*pow(dc.hx*dc.hy, 0.5);
      gsl_vector_scale (phi_r, 1.0/norm2D);
  
      gams_pqr = 0.0;
      for (n= 3; n < dc.Nfx*dc.Nfy; n++){				//excluding first 3 with zero eigen values
        gsl_matrix_get_col (psi_n, Fevec, n);
        norm2D  = gsl_blas_dnrm2(psi_n)*pow(dc.hx*dc.hy, 0.5);
        gsl_vector_scale (psi_n, 1.0/norm2D);
        xeta_n4 = gsl_vector_get(Feval, n);
         
        Hn_pq   = GalerkinProj(psi_n, phi_p, phi_q, dc);
        Hn_rs   = GalerkinProj(psi_n, phi_r, phi_s, dc);
        
        gams_pqr += Hn_pq*Hn_rs/(2.0*xeta_n4); 			// the unit of gamma here is um^-4 and needs to be converted to A^-4
      }
      gsl_matrix_set(s->alphamat, i, 0, gams_pqr*mc.Et/mc.rho*dc.gamunitconv);
      gsl_matrix_set(s->alphamat, i, 1, gams_pqr*pow(dc.Lx*dc.Ly, 3.0));
      printf("#          %03d(%02dx%02d)\t%03d(%02dx%02d)\t%03d(%02dx%02d)\t%03d(%02dx%02d)\t%f\t%5.5e\tthread:%02d\titer:%d\n", 
		  sind, mords, nords, pind, mordp, nordp, qind, mordq, nordq, rind, mordr, nordr,
		  gams_pqr*pow(dc.Lx*dc.Ly, 3.0), gams_pqr*mc.Et/mc.rho*dc.gamunitconv, omp_get_thread_num( ), i);

      // searching for same ids
      for (j = i+1; j < s->IRcou; j++){
        if ((int) gsl_matrix_get(s->IRs_pqrcombmatsorted, j, 4)==idat){
          sind = (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, j, 0);
          pind = (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, j, 1);
          qind = (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, j, 2);
          rind = (int) gsl_matrix_get(s->IRs_pqrcombmatsorted, j, 3);
          gsl_matrix_set(s->alphamat, j, 0, gams_pqr*mc.Et/mc.rho*dc.gamunitconv);
	  gsl_matrix_set(s->alphamat, j, 1, gams_pqr*pow(dc.Lx*dc.Ly, 3.0));
	  printf("#          %03d(%02dx%02d)\t%03d(%02dx%02d)\t%03d(%02dx%02d)\t%03d(%02dx%02d)\t%f\t%5.5e\tthread:%02d\titer:%d\n",
	              sind, mords, nords, pind, mordp, nordp, qind, mordq, nordq, rind, mordr, nordr,
	              gams_pqr*pow(dc.Lx*dc.Ly, 3.0), gams_pqr*mc.Et/mc.rho*dc.gamunitconv, omp_get_thread_num( ), i);
	}
      }
    }
    free(phi_s);
    free(phi_p);
    free(phi_q);
    free(phi_r);
    free(psi_n);
  }
  
  // File write !!
  // coupling file !!
  outfp = fopen("alphamat.dat", "w");
  for (i = 0; i < s->IRcou - 1; i++) {
    for (j = 0; j < 2 - 1; j++) {
      fprintf(outfp, "%5.5e\t", gsl_matrix_get(s->alphamat, i, j));
    }
    fprintf(outfp, "%5.5e\n", gsl_matrix_get(s->alphamat, i, j));
  }
  for (j = 0; j < 2 - 1; j++) {
    fprintf(outfp, "%5.5e\t", gsl_matrix_get(s->alphamat, i, j));
  }
  fprintf(outfp, "%5.5e\n", gsl_matrix_get(s->alphamat, i, j));
  fclose(outfp);
 
  return;
}
