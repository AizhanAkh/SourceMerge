/** @file nonlinear_pt.c Documented nonlinear perturbation theory module
 *
 * Mikhail M. Ivanov, 17.06.2018
 *
 * Here is a module doing one-loop perturbation theory calculations for dark matter and biased tracers
 * with IR resummation in real space. We use the FFTLog for a fast evaluation of PT loop integrals in the EdS approximation.
 * TODO: exact time dependence, redshift space, one-loop bispecrum and two-loop power spectrum.
 *
 */
#include "nonlinear_pt.h"


int nonlinear_pt_init(
                   struct precision *ppr,
                   struct background *pba,
                   struct thermo *pth,
                   struct perturbs *ppt,
                   struct primordial *ppm,
                   struct nonlinear_pt *pnlpt
                   ) {

  int index_ncdm;
  int index_k;
  int index_tau;
  double *pk_l;
  double *pk_nl;
  double *pk_Id2d2;
  double *pk_Id2;
  double *pk_IG2;
  double *pk_Id2G2;
  double *pk_IG2G2;
  double *pk_IFG2;
  double *pk_CTR;
  double *lnk_l;
  double *lnpk_l;
  double *ddlnpk_l;
  double *tau_req;
  double *deltatau;
  short print_warning=_FALSE_;
  double * pvecback;
  int last_index;
  double a,z;

  /** Summary
   *
   * (a) First deal with the case where non non-linear corrections requested */

  if (pnlpt->method == nlpt_none) {
    if (pnlpt->nonlinear_pt_verbose > 0)
      printf("No non-linear spectra requested. Nonlinear PT module skipped.\n");
  }

  /** (b) Compute the non-linear spectrum */

  else if (pnlpt->method == nlpt_spt) {
    if (pnlpt->nonlinear_pt_verbose > 0)
//      printf("Computing non-linear matter and bias power spectra in perturbation theory.\n");


      
    /** - copy list of (k,tau) from perturbation module */

    pnlpt->k_size = ppt->k_size[ppt->index_md_scalars];
    class_alloc(pnlpt->k,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    for (index_k=0; index_k<pnlpt->k_size; index_k++)
      pnlpt->k[index_k] = ppt->k[ppt->index_md_scalars][index_k];

      pnlpt->tau_size = ppt->tau_size;
      class_alloc(pnlpt->tau,pnlpt->tau_size*sizeof(double),pnlpt->error_message);
      for (index_tau=0; index_tau<pnlpt->tau_size; index_tau++){
          pnlpt->tau[index_tau] = ppt->tau_sampling[index_tau];
      }
       
      
      // Here's a new tau array (TODO: change the way the linear power spectrum is called)
      
      /*
      pnlpt->tau_size = pnlpt->z_pk_num;
      
    class_alloc(pnlpt->tau,pnlpt->tau_size*sizeof(double),pnlpt->error_message);
      
      for (index_tau=0; index_tau<pnlpt->tau_size; index_tau++){
          
          class_call(background_tau_of_z(pba,pnlpt->z_pk[pnlpt->tau_size - index_tau - 1],&pnlpt->tau[index_tau]),
                     pba->error_message,
                     pnlpt->error_message);
          printf("%lf %lf\n",pnlpt->tau[index_tau],pnlpt->z_pk[pnlpt->tau_size - index_tau - 1]);
          
      }
       */
    
//    class_alloc(pnlpt->pk_nl_out,pnlpt->tau_size*pnlpt->k_size*sizeof(double),pnlpt->error_message);
      
    class_alloc(pnlpt->nl_corr_density,pnlpt->tau_size*pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pnlpt->nl_corr_Id2d2,pnlpt->tau_size*pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pnlpt->nl_corr_Id2,pnlpt->tau_size*pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pnlpt->nl_corr_IG2,pnlpt->tau_size*pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pnlpt->nl_corr_Id2G2,pnlpt->tau_size*pnlpt->k_size*sizeof(double),pnlpt->error_message);
       class_alloc(pnlpt->nl_corr_IG2G2,pnlpt->tau_size*pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pnlpt->nl_corr_IFG2,pnlpt->tau_size*pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pnlpt->nl_corr_CTR,pnlpt->tau_size*pnlpt->k_size*sizeof(double),pnlpt->error_message);
      
//    class_alloc(pnlpt->k_nl,pnlpt->tau_size*sizeof(double),pnlpt->error_message);

    class_alloc(pk_l,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_nl,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      
    class_alloc(pk_Id2d2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_Id2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_IG2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_Id2G2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_IG2G2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_IFG2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_CTR,pnlpt->k_size*sizeof(double),pnlpt->error_message);

    class_alloc(lnk_l,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(lnpk_l,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(ddlnpk_l,pnlpt->k_size*sizeof(double),pnlpt->error_message);

      
    //  double tau_req[pnlpt->z_pk_num];
    //  double deltatau[pnlpt->z_pk_num];
      

      class_alloc(tau_req,pnlpt->z_pk_num * sizeof(double),pnlpt->error_message);
      class_alloc(deltatau,pnlpt->z_pk_num * sizeof(double),pnlpt->error_message);

      
      for (size_t i=0; i<pnlpt->z_pk_num; i++) {
          //                        printf("redshift requested %lf\n", pnlpt->z_pk[i]);
          class_call(background_tau_of_z(pba,pnlpt->z_pk[i],&tau_req[i]),
                     pba->error_message,
                     pnlpt->error_message);
      }
      
      /** Inputing the PT matrices */
      
      class_alloc(pnlpt->M22_oneline,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)) * sizeof(double),pnlpt->error_message);
      
      FILE *myFile22;
      
      myFile22 = fopen("pt_matrices/M22oneline.dat", "r");
      
      if (myFile22 == NULL){
          printf("Error Reading File M22oneline.dat\n");
          exit (0);
      }
      
      for (size_t i=0; i < (ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2; i++){
          fscanf(myFile22, "%lf", &pnlpt->M22_oneline[i]);
      }
      
      for (size_t i=(ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2; i < (ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2); i++){
          fscanf(myFile22, "%lf", &pnlpt->M22_oneline[i]);
      }
      fclose(myFile22);
      
      
      class_alloc(pnlpt->M13_oneline,((ppr->nmax_nlpt+1)*2) * sizeof(double),pnlpt->error_message);
      
      FILE *myFile13;
      
      myFile13 = fopen("pt_matrices/M13oneline.dat", "r");
      
      if (myFile13 == NULL){
          printf("Error Reading File M13oneline.dat\n");
          exit (0);
      }
      
      for (size_t i=0; i < (ppr->nmax_nlpt+1); i++){
          fscanf(myFile13, "%lf", &pnlpt->M13_oneline[i]);
      }
      
      for (size_t i=(ppr->nmax_nlpt+1); i < 2*(ppr->nmax_nlpt+1); i++){
          fscanf(myFile13, "%lf", &pnlpt->M13_oneline[i]);
      }
      fclose(myFile13);
      
      
      
      class_alloc(pnlpt->M22basic_oneline,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)) * sizeof(double),pnlpt->error_message);
      
      FILE *myFile22_basic;
      
      myFile22_basic = fopen("pt_matrices/M22basiconeline.dat", "r");
      
      if (myFile22_basic == NULL){
          printf("Error Reading File M22basiconeline.dat\n");
          exit (0);
      }
      
      for (size_t i=0; i < (ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2; i++){
          fscanf(myFile22_basic, "%lf", &pnlpt->M22basic_oneline[i]);
      }
      
      for (size_t i=(ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2; i < (ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2); i++){
          fscanf(myFile22_basic, "%lf", &pnlpt->M22basic_oneline[i]);
      }
      fclose(myFile22_basic);
      
      
      
      class_alloc(pnlpt->IFG2_oneline,((ppr->nmax_nlpt+1)*2) * sizeof(double),pnlpt->error_message);
      
      FILE *myFile_IFG2;
      
//      myFile_IFG2 = fopen("/Users/michalychforever/Dropbox/Marko/IFG2oneline.dat", "r");
      myFile_IFG2 = fopen("pt_matrices/IFG2oneline.dat", "r");
      
      if (myFile_IFG2 == NULL){
          printf("Error Reading File IFG2oneline.dat\n");
          exit (0);
      }
      
      for (size_t i=0; i < (ppr->nmax_nlpt+1); i++){
          fscanf(myFile_IFG2, "%lf", &pnlpt->IFG2_oneline[i]);
      }
      
      for (size_t i=(ppr->nmax_nlpt+1); i < 2*(ppr->nmax_nlpt+1); i++){
          fscanf(myFile_IFG2, "%lf", &pnlpt->IFG2_oneline[i]);
      }
      fclose(myFile_IFG2);

      


              /** - loop over time */
      
      int ptcounter = 0;
      double deltatautreshold = 0.00160733;
      double epsilon_for_logs2 = pow(10.,-6.);
      
          for (index_tau = pnlpt->tau_size-1; index_tau>=0; index_tau--) {

     /* Use this for debugging for a finite number of delta */
 //     for (index_tau = pnlpt->tau_size-1; index_tau > pnlpt->tau_size-2; index_tau--) {
          
  //        printf("pnlpt->tau[index_tau]=%le\n",pnlpt->tau[index_tau]);
        
       /*
        * In order to speed up calculation and don't bother changing the structure of other modules,
        * here I compute the 1-loop PS only for the redshifts close to the ones I'm interested in, $(z_my - z_in_the_loop)< treshold$,
        * and store the linear power spectrum otherwise. Then this table of P_NL(tau,k) is passed to the module spectra.c for 
        * interpolation over time and further output.
        */
              
              ptcounter = 0;
          
              for(size_t i=0; i<pnlpt->z_pk_num; i++){
                  deltatau[i] = fabs(log10(tau_req[i]) - log10(pnlpt->tau[index_tau]));
                  if(deltatau[i] < deltatautreshold){
                      ptcounter += 1;
                  }
                  else {
                      ptcounter += 0;
                  }
              }

      /* get P_L(k) at this time */
          
          class_call(nonlinear_pt_pk_l(ppt,ppm,pnlpt,index_tau,pk_l,lnk_l,lnpk_l,ddlnpk_l),
                     pnlpt->error_message,
                     pnlpt->error_message);

       /* get P_NL(k) at this time */
      if (print_warning == _FALSE_ && ptcounter > 0) {
          class_call(nonlinear_pt_loop(ppr,
                                     pba,
                                     ppm,
                                     pnlpt,
                                     pnlpt->tau[index_tau],
                                     pk_l,
                                     pk_nl,
                                     pk_Id2d2,
                                     pk_Id2,
                                     pk_IG2,
                                     pk_Id2G2,
                                     pk_IG2G2,
                                     pk_IFG2,
                                     pk_CTR,
                                     lnk_l,
                                     lnpk_l,
                                     ddlnpk_l
                                  ),
                   pnlpt->error_message,
                   pnlpt->error_message);
          
            for (index_k=0; index_k<pnlpt->k_size; index_k++) {
                pnlpt->nl_corr_density[index_tau * pnlpt->k_size + index_k] = sqrt(pk_nl[index_k]/pk_l[index_k]);
                pnlpt->nl_corr_Id2d2[index_tau * pnlpt->k_size + index_k] = sqrt(pk_Id2d2[index_k]/pk_l[index_k]);
                pnlpt->nl_corr_Id2[index_tau * pnlpt->k_size + index_k] = sqrt(pk_Id2[index_k]/pk_l[index_k]);
                pnlpt->nl_corr_IG2[index_tau * pnlpt->k_size + index_k] = sqrt(pk_IG2[index_k]/pk_l[index_k]);
                pnlpt->nl_corr_Id2G2[index_tau * pnlpt->k_size + index_k] = sqrt(pk_Id2G2[index_k]/pk_l[index_k]);
                pnlpt->nl_corr_IG2G2[index_tau * pnlpt->k_size + index_k] = sqrt(pk_IG2G2[index_k]/pk_l[index_k]);
                pnlpt->nl_corr_IFG2[index_tau * pnlpt->k_size + index_k] = sqrt(pk_IFG2[index_k]/pk_l[index_k]);
                pnlpt->nl_corr_CTR[index_tau * pnlpt->k_size + index_k] = sqrt(pk_CTR[index_k]/pk_l[index_k]);
            }


      }
          
      else {
        for (index_k=0; index_k<pnlpt->k_size; index_k++) {
          pnlpt->nl_corr_density[index_tau * pnlpt->k_size + index_k] = 1.;
          pnlpt->nl_corr_Id2[index_tau * pnlpt->k_size + index_k] = 1.;
          pnlpt->nl_corr_Id2d2[index_tau * pnlpt->k_size + index_k] = 1.;
            pnlpt->nl_corr_IG2[index_tau * pnlpt->k_size + index_k] = 1.;
            pnlpt->nl_corr_Id2G2[index_tau * pnlpt->k_size + index_k] = 1.;
            pnlpt->nl_corr_IG2G2[index_tau * pnlpt->k_size + index_k] = 1.;
            pnlpt->nl_corr_IFG2[index_tau * pnlpt->k_size + index_k] = 1.;
            pnlpt->nl_corr_CTR[index_tau * pnlpt->k_size + index_k] = 1.;
        }
      }
    }
      
    free(pk_l);
    free(pk_nl);
    free(pk_Id2d2);
    free(pk_IG2);
    free(pk_Id2G2);
    free(pk_IG2G2);
    free(pk_IFG2);
    free(pk_Id2);
    free(pk_CTR);
    free(lnk_l);
    free(lnpk_l);
    free(ddlnpk_l);
    free(deltatau);
    free(tau_req);
      
//      printf(" 'nonlinear_pt_init' module executed successfully\n");

  }
  
  else {
    class_stop(pnlpt->error_message,
               "Your non-linear method variable is set to %d, out of the range defined in nonlinear_pt.h",pnlpt->method);
  }
    
    //end of the loop over time before 'else' above

  return _SUCCESS_;
    
}


//end of non-linear init

int nonlinear_pt_free(
                   struct nonlinear_pt *pnlpt
                   ) {

  if (pnlpt->method > nlpt_none) {

    if (pnlpt->method == nlpt_spt) {
      /* free here */
        
       free(pnlpt->M13_oneline);
        free(pnlpt->M22_oneline);
        free(pnlpt->M22basic_oneline);
        free(pnlpt->IFG2_oneline);
      free(pnlpt->k);
      free(pnlpt->tau);
      free(pnlpt->nl_corr_density);
      free(pnlpt->nl_corr_Id2d2);
     free(pnlpt->nl_corr_Id2);
        free(pnlpt->nl_corr_IG2);
        free(pnlpt->nl_corr_Id2G2);
        free(pnlpt->nl_corr_IG2G2);
        free(pnlpt->nl_corr_IFG2);
        free(pnlpt->nl_corr_CTR);
 //       printf(" 'nonlinear_pt_free' module executed successfully\n");

    }
  }

  return _SUCCESS_;

}

// This function executes the linear power spectrum
int nonlinear_pt_pk_l(
                   struct perturbs *ppt,
                   struct primordial *ppm,
                   struct nonlinear_pt *pnlpt,
                   int index_tau,
                   double *pk_l,
                   double *lnk,
                   double *lnpk,
                   double *ddlnpk) {

  int index_md;
  int index_k;
  int index_ic1,index_ic2,index_ic1_ic2;
  double * primordial_pk;
  double source_ic1,source_ic2;

  index_md = ppt->index_md_scalars;

  class_alloc(primordial_pk,ppm->ic_ic_size[index_md]*sizeof(double),pnlpt->error_message);

  for (index_k=0; index_k<pnlpt->k_size; index_k++) {

    class_call(primordial_spectrum_at_k(ppm,
                                        index_md,
                                        linear,
                                        pnlpt->k[index_k],
                                        primordial_pk),
               ppm->error_message,
               pnlpt->error_message);

    pk_l[index_k] = 0;

    /* part diagonal in initial conditions */
    for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {

      index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_md]);

      source_ic1 = ppt->sources[index_md]
        [index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
        [index_tau * ppt->k_size[index_md] + index_k];
        
       // source_ic1 are transfer functions
        
      pk_l[index_k] += 2.*_PI_*_PI_/pow(pnlpt->k[index_k],3)
        *source_ic1*source_ic1
        *primordial_pk[index_ic1_ic2];
    }

      
    /* part non-diagonal in initial conditions */
    for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
      for (index_ic2 = index_ic1+1; index_ic2 < ppm->ic_size[index_md]; index_ic2++) {

        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,ppm->ic_size[index_md]);

        if (ppm->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

          source_ic1 = ppt->sources[index_md]
            [index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
            [index_tau * ppt->k_size[index_md] + index_k];

          source_ic2 = ppt->sources[index_md]
            [index_ic2 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
            [index_tau * ppt->k_size[index_md] + index_k];

          pk_l[index_k] += 2.*2.*_PI_*_PI_/pow(pnlpt->k[index_k],3)
            *source_ic1*source_ic2
            *primordial_pk[index_ic1_ic2]; // extra 2 factor (to include the symmetric term ic2,ic1)

        }
      }
    }

    lnk[index_k] = log(pnlpt->k[index_k]);
    lnpk[index_k] = log(pk_l[index_k]);
 //   printf("%e  %e\n",pnlpt->k[index_k],pk_l[index_k]);
  }

    class_call(array_spline_table_columns(lnk,
                                          pnlpt->k_size,
                                          lnpk,
                                          1,
                                          ddlnpk,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

  free(primordial_pk);

  return _SUCCESS_;

} // end of non-linear_Pk_l. This function executes the linear power spectrum


/* beginning of my own function */
 
 int nonlinear_pt_loop(
                       struct precision *ppr,
                       struct background *pba,
                       struct primordial *ppm,
                       struct nonlinear_pt *pnlpt,
                       double tau,
                       double *pk_l,
                       double *pk_nl,
                       double *pk_Id2d2,
                       double *pk_Id2,
                       double *pk_IG2,
                       double *pk_Id2G2,
                       double *pk_IG2G2,
                       double *pk_IFG2,
                       double *pk_CTR,
                       double *lnk_l,
                       double *lnpk_l,
                       double *ddlnpk_l
                       ) {

int index_k;
double pk_lin,pk_quasi,pk_halo,rk;
double sigmav;
     
double * pvecback;

int last_index;
double sum0;
double anorm;

double *integrand_array;
int integrand_size;
int index_ia_k;
int index_ia_pk;
int index_ia_sum;
int index_ia_ddsum;

int ia_size;
int index_ia;

double k_integrand;
double lnpk_integrand;
     
//* halofit_found_k_max = ok;
     
class_alloc(pvecback,pba->bg_size*sizeof(double),pnlpt->error_message);
     
anorm    = 1./(2*pow(_PI_,2));
 

index_ia=0;
class_define_index(index_ia_k,     _TRUE_,index_ia,1);
class_define_index(index_ia_pk,    _TRUE_,index_ia,1);
class_define_index(index_ia_sum,   _TRUE_,index_ia,1);
class_define_index(index_ia_ddsum, _TRUE_,index_ia,1);
ia_size = index_ia;

//M: here the integrand_array is tabulated. First one computes the size of the integrand (number of points I presume)

integrand_size=(int)(log(pnlpt->k[pnlpt->k_size-1]/pnlpt->k[0])/log(10.)*ppr->halofit_k_per_decade)+1;
class_alloc(integrand_array,integrand_size*ia_size*sizeof(double),pnlpt->error_message);

/* we fill integrand_array with values of k and P(k) using interpolation */

last_index=0;

for (index_k=0; index_k < integrand_size; index_k++) {
    
    k_integrand=pnlpt->k[0]*pow(10.,index_k/ppr->halofit_k_per_decade);
    
    class_call(array_interpolate_spline(lnk_l,
                                        pnlpt->k_size,
                                        lnpk_l,           /*what is being interpolated*/
                                        ddlnpk_l,
                                        1,
                                        log(k_integrand),
                                        &last_index,
                                        &lnpk_integrand,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    integrand_array[index_k*ia_size + index_ia_k] = k_integrand;
    integrand_array[index_k*ia_size + index_ia_pk] = exp(lnpk_integrand);
    
}

// class_call only calls a function!

class_call(background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback),
                pba->error_message,
                pnlpt->error_message);
     
printf("Computing one-loop non-linear matter and bias power spectra at z=%e\n",pba->a_today/pvecback[pba->index_bg_a]-1.);

class_call(nonlinear_pt_p13UV_integrate(
                                     pnlpt,
                                     pba,
                                     integrand_array,
                                     integrand_size,
                                     ia_size,
                                     index_ia_k,
                                     index_ia_pk,
                                     index_ia_sum,
                                     index_ia_ddsum,
                                     &sum0
                                     ),
           pnlpt->error_message,
           pnlpt->error_message);
//printf("Sigma_v=%f\n",sum0);
sigmav = sum0;
     
     
     free(pvecback);
     free(integrand_array);
     
     int Nmax = ppr->nmax_nlpt;
     double Nmaxd = Nmax * 1.;
//     double kmin = 0.0001 * pba->h;
//     double kmax = 100. * pba->h;
     
     double kmin = pnlpt->k[0];
     double kmax = pnlpt->k[pnlpt->k_size-1];
     
     
     double *js;
     class_alloc(js,(Nmax+1)*sizeof(double),pnlpt->error_message);
     double *kdisc;
     class_alloc(kdisc,Nmax * sizeof(double),pnlpt->error_message);
     double *Pbin;
     class_alloc(Pbin,Nmax * sizeof(double),pnlpt->error_message);
     double *Pdisc;
     class_alloc(Pdisc,Nmax * sizeof(double),pnlpt->error_message);
     double *Ptree;
     class_alloc(Ptree,Nmax * sizeof(double),pnlpt->error_message);
     
     
     // double kmin = pnl->k[0];
     double Delta = log(kmax / kmin) / (Nmax - 1);
     
     /*
     double js[Nmax+1], kdisc[Nmax], Pdisc[Nmax], Pbin[Nmax];
     double complex etam[Nmax+1];
      */
     
     // Here I interpolate the linear power spectrum
     
     gsl_interp_accel *acc_lnpl
     = gsl_interp_accel_alloc();
     gsl_spline *spline_lnpl
     = gsl_spline_alloc (gsl_interp_cspline, pnlpt->k_size);
     gsl_spline_init(spline_lnpl, lnk_l, lnpk_l, pnlpt->k_size);
     
     for (size_t i=0; i< Nmax; i++){
         kdisc[i] = kmin * exp(i * Delta);
         Pdisc[i] = exp(gsl_spline_eval(spline_lnpl, log(kdisc[i]), acc_lnpl));
     }
     
     if (pnlpt->irres == irres_yes) {
         
    // Here I do IR resummation
         
     printf("IR resummation performed.\n");
         
     //   IR-1) Computing the DFST-II of log(kP)
     
     int Nir = 65536;
     int Nirover2 = 32768;
     double Nird = Nir * 1.;
     double kmin2 = 0.0001 * pba->h;
     double kmax2 = 10.* pba->h;
     double *logkPdiscr;
     class_alloc(logkPdiscr,Nir * sizeof(double),pnlpt->error_message);
     double Pbin2;
     double kbin2;
         
     
     for (size_t i=0; i< Nir; i++){
         kbin2 = kmin2 + i * (kmax2 - kmin2)/(Nird-1.);
         // exp(gsl_spline_eval(spline_lnpl, log(kdisc[i]), acc_lnpl));
         Pbin2 = exp(gsl_spline_eval(spline_lnpl, log(kbin2), acc_lnpl));
         logkPdiscr[i] = log(kbin2 * Pbin2);
    //     printf("%le\n",logkPdiscr[i]);
     }
     
     gsl_spline_free (spline_lnpl);
     gsl_interp_accel_free (acc_lnpl);
     
         
     fftw_plan p_ir;
//     double out_ir[Nir];
         
         double *out_ir;
         class_alloc(out_ir,Nir * sizeof(double),pnlpt->error_message);
         
     p_ir = fftw_plan_r2r_1d(Nir, logkPdiscr, out_ir, FFTW_RODFT10, FFTW_ESTIMATE);
     fftw_execute(p_ir);
     
     
         double *cmodd;
         class_alloc(cmodd,Nirover2 * sizeof(double),pnlpt->error_message);
         double *cmeven;
         class_alloc(cmeven,Nirover2 * sizeof(double),pnlpt->error_message);
         int *ivar;
         class_alloc(ivar,Nirover2 * sizeof(int),pnlpt->error_message);
         
     for (size_t i=0; i<Nirover2; i++){
         ivar[i] = i + 1;
         cmodd[i] = out_ir[2*i];
         cmeven[i] = out_ir[2*i+1];
     }
     
     fftw_destroy_plan(p_ir);
     
     //   IR-2) Erasing the BAO bump from the odd and even Fourier harmonics and smoothy interpolating the erased Fourier coefficients
         
     int Nleft = 120;
     int Nright = 220;
     int Nthrow = Nright - Nleft;
     int Nnew = Nirover2 - Nthrow;
         double *cmoddnw;
         class_alloc(cmoddnw,Nnew * sizeof(double),pnlpt->error_message);
         double *cmevennw;
         class_alloc(cmevennw,Nnew * sizeof(double),pnlpt->error_message);
         double *inew;
         class_alloc(inew,Nnew * sizeof(double),pnlpt->error_message);
     
     for (size_t i=0; i<Nleft; i++){
         cmoddnw[i] = cmodd[i];
         cmevennw[i] = cmeven[i];
         inew[i] = i+1.;
     }
     
     for (size_t i= Nleft; i<Nnew; i++){
         cmoddnw[i] = cmodd[i + Nthrow];
         cmevennw[i] = cmeven[i + Nthrow];
         inew[i] = i+1. + Nthrow;
     }
         
     
     gsl_interp_accel *acc_cmoddnw
     = gsl_interp_accel_alloc();
     gsl_spline *spline_cmoddnw
     = gsl_spline_alloc (gsl_interp_cspline, Nnew);
     gsl_spline_init(spline_cmoddnw, inew, cmoddnw, Nnew);
     
     gsl_interp_accel *acc_cmevennw
     = gsl_interp_accel_alloc();
     gsl_spline *spline_cmevennw
     = gsl_spline_alloc (gsl_interp_cspline, Nnew);
     gsl_spline_init(spline_cmevennw, inew, cmevennw, Nnew);
     
         double *cmnew;
         class_alloc(cmnew,Nir * sizeof(double),pnlpt->error_message);
         
     for (size_t i=0; i<Nirover2; i++ ){
         cmnew[i*2] = gsl_spline_eval(spline_cmoddnw, 1.*i + 1., acc_cmoddnw);
         cmnew[i*2+1] = gsl_spline_eval(spline_cmevennw, 1.*i + 1., acc_cmevennw);
     }
     
     //   IR-3) Inverse DST-II (= DST-III/(2*N)) and interpolating P_nw
         
     
     fftw_plan p_2;
         
         double *out_2;
         class_alloc(out_2,Nir * sizeof(double),pnlpt->error_message);
         
     p_2 = fftw_plan_r2r_1d(Nir, cmnew, out_2, FFTW_RODFT01, FFTW_ESTIMATE);
     fftw_execute(p_2);
         
         double *Pnw_ir;
         class_alloc(Pnw_ir,Nir * sizeof(double),pnlpt->error_message);
         double *knw_ir;
         class_alloc(knw_ir,Nir * sizeof(double),pnlpt->error_message);
         
     for (size_t i=0; i<Nir; i++){
         knw_ir[i] = kmin2 + i * (kmax2 - kmin2)/(Nird-1.);
         Pnw_ir[i] = exp(out_2[i]/(2*Nird))/(knw_ir[i]);
         //       printf("%le %le\n",knw[i],Pnw[i]);
     }
     
     fftw_destroy_plan(p_2);
     gsl_spline_free (spline_cmevennw);
     gsl_spline_free (spline_cmoddnw);
     gsl_interp_accel_free (acc_cmevennw);
     gsl_interp_accel_free (acc_cmoddnw);
         
     gsl_interp_accel *acc_nw
     = gsl_interp_accel_alloc();
     gsl_spline *spline_nw
     = gsl_spline_alloc (gsl_interp_cspline, Nir);
     gsl_spline_init(spline_nw, knw_ir, Pnw_ir, Nir);
         
     
     //   IR-4) Computing the BAO damping factor
     
     double SigmaBAO;
     double rbao = 110./pba->h;
     int Nint2 = 500;
//     double qint2[Nint2+1], IntegrandBAO[Nint2+1];
         
         double *qint2;
         class_alloc(qint2,Nint2 * sizeof(double),pnlpt->error_message);
         
         double *IntegrandBAO;
         class_alloc(IntegrandBAO,Nint2 * sizeof(double),pnlpt->error_message);
         
         double ks = 0.2 * pba->h;
         for (size_t i=0; i<Nint2+1; i++){
             qint2[i] = kmin2 * exp(i * log(ks/kmin2) / (Nint2));
             IntegrandBAO[i] = gsl_spline_eval(spline_nw, qint2[i], acc_nw) * (1. - 3.*sin(qint2[i] * rbao)/(qint2[i] * rbao) + 6.*(sin(qint2[i] * rbao)/pow((qint2[i] * rbao),3.) - cos(qint2[i] * rbao)/pow((qint2[i] * rbao),2.)));
         }
     
         for (size_t i=0; i<Nint2; i++){
             SigmaBAO += (log( qint2[i+1])-log( qint2[i])) *(qint2[i+1]*IntegrandBAO[i+1] + qint2[i]*IntegrandBAO[i])/2 / (6 * pow(M_PI,2));
         }
//     double SigmaBAOh = SigmaBAO * pow(pba->h,2.);
//     printf("Sigma_BAO(ks=0.2 h/Mpc)=%lf (Mpc/h)^2\n",SigmaBAOh);
     
     //   IR-5) Computing the LO IR resummed power spectrum
     
         
     
         double *Pnw;
         class_alloc(Pnw,Nmax * sizeof(double),pnlpt->error_message);
         
         double *Pw;
         class_alloc(Pw,Nmax * sizeof(double),pnlpt->error_message);
         
     
         for (size_t i=0; i<Nmax; i++){
         
             if (kdisc[i]<= kmax2 && kdisc[i] >= kmin2) {
                 Pnw[i] = gsl_spline_eval(spline_nw, kdisc[i], acc_nw);
                 Pw[i] = Pdisc[i] - Pnw[i];
                 Pbin[i] = Pnw[i] + Pw[i] * exp(-SigmaBAO * pow(kdisc[i],2.));
             }
         
         else {
                Pnw[i] = Pdisc[i];
                Pw[i] = 0.;
                Pbin[i] = Pdisc[i];
         }
         Ptree[i] = Pnw[i] + Pw[i] * exp(-SigmaBAO * pow(kdisc[i],2.))*(1. + SigmaBAO * pow(kdisc[i],2.));
     }
     
     gsl_spline_free (spline_nw);
     gsl_interp_accel_free (acc_nw);
         
         
         free(out_ir);
         free(out_2);
         free(cmnew);
         free(logkPdiscr);
         free(cmeven);
         free(cmodd);
         free(ivar);
         free(inew);
         free(cmoddnw);
         free(cmevennw);
         free(knw_ir);
         free(Pnw_ir);
         free(qint2);
         free(IntegrandBAO);
         free(Pw);
         free(Pnw);
         
     
} /* End of IR resummation conditional expression */
     
     else{
         
         gsl_spline_free (spline_lnpl);
         gsl_interp_accel_free (acc_lnpl);
         
         printf("IR resummation skipped.\n");
         for (size_t i=0; i<Nmax; i++){
          Pbin[i] = Pdisc[i];
          Ptree[i] = Pbin[i];
        }
       }
     
     // here we compute the FFT coefficients
     
     double complex *etam;
     class_alloc(etam,(Nmax+1)*sizeof(complex double),pnlpt->error_message);
     
     
     double b = -0.3;
     for (size_t i=0; i< Nmax +1 ; i++){
         js[i] = i - Nmaxd/2;
         etam[i] = b + 2 * M_PI * _Complex_I * js[i]/Nmaxd / Delta ;
     }
     
      fftw_complex *in, *out;
      fftw_plan p;
     
     in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);
     out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);
     
     for (size_t i=0; i< Nmax ; i++){
         in[i][0] = Pbin[i]* exp(-1.* i * b* Delta);
         in[i][1] = 0;
     }
     
      p = fftw_plan_dft_1d(Nmax, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
      fftw_execute(p);
     
     
     /*
      printf("Start of cn's\n");
      
      for (size_t i=0; i< Nmax+1; i++){
      printf("%f %f\n",out[i][0],out[i][1]);
      }
      
      printf("End of cn's\n");
      */
     
     
 //    double complex *cmsym = malloc((Nmax+1) * sizeof(complex double));
     
     double complex *cmsym;
     class_alloc(cmsym,(Nmax+1)*sizeof(complex double),pnlpt->error_message);
     
//     double complex cmsym[Nmax+1];
     for (size_t i=0; i< Nmax+1; i++){
         if (i < Nmax/2) {
            cmsym[i]= cpow(kmin,-etam[i]) * (out[Nmax/2 - i][0] - _Complex_I * out[Nmax/2 - i][1])/Nmaxd;
         }
         else {
            cmsym[i]= cpow(kmin,-etam[i]) * (out[i - Nmax/2][0] + _Complex_I * out[i - Nmax/2][1])/Nmaxd;
         }
     }
     
     cmsym[0] = cmsym[0]/2.;
     cmsym[Nmax] = cmsym[Nmax]/2.;
     
     
     fftw_free(in);
     fftw_free(out);
     fftw_destroy_plan(p);
     
     
     
     
     /*
     printf("Start of cnsym's\n");
     
     for (size_t i=0; i< Nmax+1; i++){
         printf("%le  %le\n",creal(cmsym[i]),cimag(cmsym[i]));
     }
     
     printf("End of cnsym's\n");
      */
     
     
     
// here we input the precomputed PT matrices
     
     double cutoff = 3. * pba->h;
//     double P13[Nmax], P13UV[Nmax], P1loop[Nmax];
     
     double *P13;
     double *P13UV;
     double *P1loop;
     double complex *f13;
     class_alloc(P13,Nmax*sizeof(double),pnlpt->error_message);
     class_alloc(P13UV,Nmax*sizeof(double),pnlpt->error_message);
     class_alloc(P1loop,Nmax*sizeof(double),pnlpt->error_message);
     class_alloc(f13,Nmax*sizeof(complex double),pnlpt->error_message);
     
     double complex *f22;
     double *P22;
     double *P_CTR;
     class_alloc(f22,Nmax*sizeof(complex double),pnlpt->error_message);
     class_alloc(P22,Nmax*sizeof(double),pnlpt->error_message);
     class_alloc(P_CTR,Nmax*sizeof(double),pnlpt->error_message);
     
 //    double P22[Nmax];
     int count=0;
     for (size_t j=0; j < Nmax; j++){
         f22[j] = 0.;
         
         count=0;
         
         for (size_t i=0; i < Nmax+1; i++){
             
             for (size_t l=0; l <= i; l++){
                 
                 if (l != i){
                     f22[j] += 2. * (cmsym[l]* cpow(kdisc[j], etam[l])) * (cmsym[i]* cpow(kdisc[j], etam[i])) * (pnlpt->M22_oneline[count] + _Complex_I * pnlpt->M22_oneline[count + (Nmax+1)*(Nmax+2)/2]);
                                      }
                 else {
                     f22[j] += cpow((cmsym[i]* cpow(kdisc[j], etam[i])),2.) * (pnlpt->M22_oneline[count] + _Complex_I * pnlpt->M22_oneline[count + (Nmax+1)*(Nmax+2)/2]);
                 }
                 count++;
             }
         }
         P22[j] = creal(cpow(kdisc[j], 3) * f22[j] * exp(-pow(kdisc[j]/cutoff, 6)) );
//         printf("%le %le\n",kdisc[j],P22[j]);
     }
     
     
     for (size_t j=0; j < Nmax; j++){
         f13[j] = 0.;
         for (size_t i=0; i < Nmax+1; i++){
             f13[j] += cmsym[i]* cpow(kdisc[j], etam[i]) * (pnlpt->M13_oneline[i]+ _Complex_I * pnlpt->M13_oneline[i + Nmax+1]);
         }
         
         P13UV[j] = - 61. * Pbin[j] * pow(kdisc[j], 2.) * sigmav / 105.;
         P13[j] = (creal(cpow(kdisc[j], 3) * f13[j] * Pbin[j]) + P13UV[j]) * exp(-pow(kdisc[j]/cutoff, 6.));
         P1loop[j] = Ptree[j] +(P13[j] + P22[j]);
         P_CTR[j] = kdisc[j] * kdisc[j] * Pbin[j];
     }
     
     // Here I interpolate the non-linear power spectrum to pass it further
     
     gsl_interp_accel *acc
     = gsl_interp_accel_alloc();
     gsl_spline *spline
     = gsl_spline_alloc (gsl_interp_cspline, Nmax);
     gsl_spline_init(spline, kdisc, P1loop, Nmax);
     
     
     gsl_interp_accel *acc_CTR
     = gsl_interp_accel_alloc();
     gsl_spline *spline_CTR
     = gsl_spline_alloc (gsl_interp_cspline, Nmax);
     gsl_spline_init(spline_CTR, kdisc, P_CTR, Nmax);
     
     
     for (index_k=0; index_k < pnlpt->k_size; index_k++){
         pk_lin = pk_l[index_k]*pow(pnlpt->k[index_k],3)*anorm;
         
//         if (pnlpt->k[index_k]>= kmin && pnlpt->k[index_k]<= kmax){
             pk_nl[index_k] = gsl_spline_eval(spline, pnlpt->k[index_k], acc);
             pk_CTR[index_k] = gsl_spline_eval(spline_CTR, pnlpt->k[index_k], acc_CTR);
//         }
         
//         else{
//             pk_nl[index_k] = pk_l[index_k];
//         }
     }
     gsl_spline_free (spline);
     gsl_interp_accel_free (acc);
     
     gsl_spline_free (spline_CTR);
     gsl_interp_accel_free (acc_CTR);
     
     free(etam);
     free(cmsym);
     free(P22);
     free(P13);
     free(P13UV);
     free(P_CTR);
     free(P1loop);
     free(f13);
     free(f22);
 //    free(M13_oneline);
 //    free(M22_oneline);
     free(Ptree);
     
     
     /* Here I compute the power spectra for biased tracers. For this reason I have to compute the FFTLog coefficients for a new bias b2 */

     double complex *etam2;
     class_alloc(etam2,(Nmax+1)*sizeof(complex double),pnlpt->error_message);
     
     double b2 = -1.6000001;
     for (size_t i=0; i< Nmax +1 ; i++){
         etam2[i] = b2 + 2 * M_PI * _Complex_I * js[i]/Nmaxd / Delta ;
     }
     
     
     fftw_complex *in_bias, *out_bias;
     fftw_plan p_bias;
     
     in_bias = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);
     out_bias = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);
     
     for (size_t i=0; i< Nmax ; i++){
         in_bias[i][0] = Pbin[i]* exp(-1.* i * b2* Delta);
         in_bias[i][1] = 0;
     }
     
     p_bias = fftw_plan_dft_1d(Nmax, in_bias, out_bias, FFTW_FORWARD, FFTW_ESTIMATE);
     fftw_execute(p_bias);
     
     
     double complex *cmsym2;
     class_alloc(cmsym2,(Nmax+1)*sizeof(complex double),pnlpt->error_message);
     
     for (size_t i=0; i< Nmax+1; i++){
         if (i < Nmax/2) {
             cmsym2[i]= cpow(kmin,-etam2[i]) * (out_bias[Nmax/2 - i][0] - _Complex_I * out_bias[Nmax/2 - i][1])/Nmaxd;
         }
         else {
             cmsym2[i]= cpow(kmin,-etam2[i]) * (out_bias[i - Nmax/2][0] + _Complex_I * out_bias[i - Nmax/2][1])/Nmaxd;
         }
     }
     
     
     cmsym2[0] = cmsym2[0]/2.;
     cmsym2[Nmax] = cmsym2[Nmax]/2.;
     
     
     fftw_free(in_bias);
     fftw_free(out_bias);
     fftw_destroy_plan(p_bias);
     
          /* Here I compute Id2d2 */

     
     double complex *f22_Id2d2;
     double *P_Id2d2;
     class_alloc(f22_Id2d2,Nmax*sizeof(complex double),pnlpt->error_message);
     class_alloc(P_Id2d2,Nmax*sizeof(double),pnlpt->error_message);
     
     double epsilon_for_logs = 1.*pow(10.,-6.);
     
     int count2=0;
     for (size_t j=0; j < Nmax; j++){
         f22_Id2d2[j] = 0.;
         
         count2=0;
         
         for (size_t i=0; i < Nmax+1; i++){
             
             for (size_t l=0; l <= i; l++){
                 
                 if (l != i){
                     f22_Id2d2[j] += 2. * (cmsym2[l]* cpow(kdisc[j], etam2[l])) * (cmsym2[i]* cpow(kdisc[j], etam2[i])) * 2. * (pnlpt->M22basic_oneline[count2] + _Complex_I * pnlpt->M22basic_oneline[count2 + (Nmax+1)*(Nmax+2)/2]);
                 }
                 else {
                     f22_Id2d2[j] += cpow((cmsym2[i]* cpow(kdisc[j], etam2[i])),2.) * 2. * (pnlpt->M22basic_oneline[count2] + _Complex_I * pnlpt->M22basic_oneline[count2 + (Nmax+1)*(Nmax+2)/2]);
                 }
                 count2++;
             }
         }
         P_Id2d2[j] = fabs(creal(cpow(kdisc[j], 3) * f22_Id2d2[j]) - creal(cpow(kdisc[0], 3) * f22_Id2d2[0]) + epsilon_for_logs);
     }
     
      gsl_interp_accel *acc_PId2d2
     = gsl_interp_accel_alloc();
     gsl_spline *spline_PId2d2
     = gsl_spline_alloc (gsl_interp_cspline, Nmax);
     gsl_spline_init(spline_PId2d2, kdisc, P_Id2d2, Nmax);

     
     
     for (index_k=0; index_k < pnlpt->k_size; index_k++){
                      pk_Id2d2[index_k] = gsl_spline_eval(spline_PId2d2, pnlpt->k[index_k], acc_PId2d2);
     }
     
     gsl_spline_free (spline_PId2d2);
     gsl_interp_accel_free (acc_PId2d2);
     
     free(f22_Id2d2);
     free(P_Id2d2);
     
     /* Here I compute Id2 */
     
     
     double complex *f22_Id2;
     double *P_Id2;
     class_alloc(f22_Id2,Nmax*sizeof(complex double),pnlpt->error_message);
     class_alloc(P_Id2,Nmax*sizeof(double),pnlpt->error_message);

     
     count2=0;
     
     for (size_t j=0; j < Nmax; j++){
         f22_Id2[j] = 0.;
         
         count2=0;
         
         for (size_t i=0; i < Nmax+1; i++){
             
             for (size_t l=0; l <= i; l++){
                 
                 if (l != i){
                     f22_Id2[j] += 2. * (cmsym2[l]* cpow(kdisc[j], etam2[l])) * (cmsym2[i]* cpow(kdisc[j], etam2[i])) * (pnlpt->M22basic_oneline[count2] + _Complex_I * pnlpt->M22basic_oneline[count2 + (Nmax+1)*(Nmax+2)/2]) * ((3. + etam2[i] + etam2[l])*(4. + 3.5 *(etam2[i] + etam2[l]))/(14.*(-0.5*etam2[l])*(-0.5 *etam2[i])));
                 }
                 else {
                     f22_Id2[j] += cpow((cmsym2[i]* cpow(kdisc[j], etam2[i])),2.) * (pnlpt->M22basic_oneline[count2] + _Complex_I * pnlpt->M22basic_oneline[count2 + (Nmax+1)*(Nmax+2)/2])* ((3. + etam2[i] + etam2[i])*(4. + 3.5 *(etam2[i] + etam2[i]))/(14. *(-0.5*etam2[i])*(-0.5 *etam2[i])));
                 }
                 count2++;
             }
         }
         P_Id2[j] = creal(cpow(kdisc[j], 3) * f22_Id2[j]);
//                        printf("%le %le\n",kdisc[j],P_Id2[j]);
     }
     
     // Here I interpolate the non-linear power spectrum to pass it further
     
     gsl_interp_accel *acc_PId2
     = gsl_interp_accel_alloc();
     gsl_spline *spline_PId2
     = gsl_spline_alloc (gsl_interp_cspline, Nmax);
     gsl_spline_init(spline_PId2, kdisc, P_Id2, Nmax);
     
     
     for (index_k=0; index_k < pnlpt->k_size; index_k++){
         
//         if (pnlpt->k[index_k]>= kmin && pnlpt->k[index_k]<= kmax){
             pk_Id2[index_k] = pk_l[index_k] + gsl_spline_eval(spline_PId2, pnlpt->k[index_k], acc_PId2);
//         }
         
//         else{
//         pk_Id2[index_k] =  pk_l[index_k];
//         }
     }
     
     gsl_spline_free (spline_PId2);
     gsl_interp_accel_free (acc_PId2);

     
//     free(Id2_oneline);
     free(f22_Id2);
     free(P_Id2);
     
     
     /* Here I compute IG2 */
     
     double complex *f22_IG2;
     double *P_IG2;
     class_alloc(f22_IG2,Nmax*sizeof(complex double),pnlpt->error_message);
     class_alloc(P_IG2,Nmax*sizeof(double),pnlpt->error_message);
     
     for (size_t j=0; j < Nmax; j++){
         f22_IG2[j] = 0.;
         
         count2=0;
         
         for (size_t i=0; i < Nmax+1; i++){
             
             for (size_t l=0; l <= i; l++){
                 
                 if (l != i){
                     f22_IG2[j] += 2. * (cmsym2[l]* cpow(kdisc[j], etam2[l])) * (cmsym2[i]* cpow(kdisc[j], etam2[i])) * (pnlpt->M22basic_oneline[count2] + _Complex_I * pnlpt->M22basic_oneline[count2 + (Nmax+1)*(Nmax+2)/2])*(-1.*(3.+etam2[i]+etam2[l])*(1.+etam2[i]+etam2[l])*(6. - 3.5 *(etam2[i]+etam2[l]))/(28.*(1.-0.5*etam2[i])*(1.-0.5*etam2[l])*(-0.5*etam2[i])*(-0.5*etam2[l])));
                 }
                 else {
                     f22_IG2[j] += cpow((cmsym2[i]* cpow(kdisc[j], etam2[i])),2.) * (pnlpt->M22basic_oneline[count2] + _Complex_I * pnlpt->M22basic_oneline[count2 + (Nmax+1)*(Nmax+2)/2])*(-1.*(3.+etam2[i]+etam2[i])*(1.+etam2[i]+etam2[i])*(6. - 3.5 *(etam2[i]+etam2[i]))/(28.*(1.-0.5*etam2[i])*(1.-0.5*etam2[i])*(-0.5*etam2[i])*(-0.5*etam2[i])));
                 }
                 count2++;
             }
         }
         P_IG2[j] = fabs(creal(cpow(kdisc[j], 3) * f22_IG2[j]));
         //                        printf("%le %le\n",kdisc[j],P_Id2[j]);
     }
     
     
     //               printf("everything is OK here 2\n");
     
     // Here I interpolate the non-linear power spectrum to pass it further
     
     gsl_interp_accel *acc_PIG2
     = gsl_interp_accel_alloc();
     gsl_spline *spline_PIG2
     = gsl_spline_alloc (gsl_interp_cspline, Nmax);
     gsl_spline_init(spline_PIG2, kdisc, P_IG2, Nmax);
     
     for (index_k=0; index_k < pnlpt->k_size; index_k++){
        
         pk_IG2[index_k] = gsl_spline_eval(spline_PIG2, pnlpt->k[index_k], acc_PIG2);
     }
     
     gsl_spline_free (spline_PIG2);
     gsl_interp_accel_free (acc_PIG2);
     
     free(f22_IG2);
     free(P_IG2);


     
     
     /* Here I compute Id2G2 */
     
     double complex *f22_Id2G2;
     double *P_Id2G2;
     class_alloc(f22_Id2G2,Nmax*sizeof(complex double),pnlpt->error_message);
     class_alloc(P_Id2G2,Nmax*sizeof(double),pnlpt->error_message);
     
     for (size_t j=0; j < Nmax; j++){
         f22_Id2G2[j] = 0.;
         
         count2=0;
         
         for (size_t i=0; i < Nmax+1; i++){
             
             for (size_t l=0; l <= i; l++){
                 
                 if (l != i){
                     f22_Id2G2[j] += 2. * (cmsym2[l]* cpow(kdisc[j], etam2[l])) * (cmsym2[i]* cpow(kdisc[j], etam2[i])) * (pnlpt->M22basic_oneline[count2] + _Complex_I * pnlpt->M22basic_oneline[count2 + (Nmax+1)*(Nmax+2)/2])*((3.+etam2[i]+etam2[l])/((-0.5*etam2[i])*(-0.5*etam2[l])));
                 }
                 else {
                     f22_Id2G2[j] += cpow((cmsym2[i]* cpow(kdisc[j], etam2[i])),2.) * (pnlpt->M22basic_oneline[count2] + _Complex_I * pnlpt->M22basic_oneline[count2 + (Nmax+1)*(Nmax+2)/2])*((3.+etam2[i]+etam2[i])/((-0.5*etam2[i])*(-0.5*etam2[i])));
                 }
                 count2++;
             }
         }
         P_Id2G2[j] = fabs(creal(cpow(kdisc[j], 3) * f22_Id2G2[j]));
     }
     
     gsl_interp_accel *acc_PId2G2
     = gsl_interp_accel_alloc();
     gsl_spline *spline_PId2G2
     = gsl_spline_alloc (gsl_interp_cspline, Nmax);
     gsl_spline_init(spline_PId2G2, kdisc, P_Id2G2, Nmax);
     
     
     for (index_k=0; index_k < pnlpt->k_size; index_k++){
         
         pk_Id2G2[index_k] = gsl_spline_eval(spline_PId2G2, pnlpt->k[index_k], acc_PId2G2);
     }
     
     gsl_spline_free (spline_PId2G2);
     gsl_interp_accel_free (acc_PId2G2);
     
     free(f22_Id2G2);
     free(P_Id2G2);
     
     
     
     /* Here I compute IG2G2 */
     
     double complex *f22_IG2G2;
     double *P_IG2G2;
     class_alloc(f22_IG2G2,Nmax*sizeof(complex double),pnlpt->error_message);
     class_alloc(P_IG2G2,Nmax*sizeof(double),pnlpt->error_message);
     
     
     for (size_t j=0; j < Nmax; j++){
         f22_IG2G2[j] = 0.;
         
         count2=0;
         
         for (size_t i=0; i < Nmax+1; i++){
             
             for (size_t l=0; l <= i; l++){
                 
                 if (l != i){
                     f22_IG2G2[j] += 2. * (cmsym2[l]* cpow(kdisc[j], etam2[l])) * (cmsym2[i]* cpow(kdisc[j], etam2[i])) * (pnlpt->M22basic_oneline[count2] + _Complex_I * pnlpt->M22basic_oneline[count2 + (Nmax+1)*(Nmax+2)/2])*((3.+etam2[i]+etam2[l])*(1.+etam2[i]+etam2[l])/((-0.5*etam2[i])*(-0.5*etam2[l])*(1.-0.5*etam2[i])*(1.-0.5*etam2[l])));
                 }
                 else {
                     f22_IG2G2[j] += cpow((cmsym2[i]* cpow(kdisc[j], etam2[i])),2.) * (pnlpt->M22basic_oneline[count2] + _Complex_I * pnlpt->M22basic_oneline[count2 + (Nmax+1)*(Nmax+2)/2])*((3.+2.*etam2[i])*(1.+2.*etam2[i])/((-0.5*etam2[i])*(-0.5*etam2[i])*(1.-0.5*etam2[i])*(1.-0.5*etam2[i])));
                 }
                 count2++;
             }
         }
         P_IG2G2[j] = fabs(creal(cpow(kdisc[j], 3) * f22_IG2G2[j]));
     }
     
     
     gsl_interp_accel *acc_PIG2G2
     = gsl_interp_accel_alloc();
     gsl_spline *spline_PIG2G2
     = gsl_spline_alloc (gsl_interp_cspline, Nmax);
     gsl_spline_init(spline_PIG2G2, kdisc, P_IG2G2, Nmax);
     
     for (index_k=0; index_k < pnlpt->k_size; index_k++){
         
         pk_IG2G2[index_k] = gsl_spline_eval(spline_PIG2G2, pnlpt->k[index_k], acc_PIG2G2);
     }
     
     gsl_spline_free (spline_PIG2G2);
     gsl_interp_accel_free (acc_PIG2G2);
     
     free(f22_IG2G2);
     free(P_IG2G2);
 
     /* Here I compute IFG2 */
     
     double complex *f13_IFG2;
     double *P_IFG2;
     class_alloc(f13_IFG2,Nmax*sizeof(complex double),pnlpt->error_message);
     class_alloc(P_IFG2,Nmax*sizeof(double),pnlpt->error_message);
     
     for (size_t j=0; j < Nmax; j++){
         f13_IFG2[j] = 0.;
         for (size_t i=0; i < Nmax+1; i++){
             f13_IFG2[j] += cmsym2[i]* cpow(kdisc[j], etam2[i]) * (pnlpt->IFG2_oneline[i]+ _Complex_I * pnlpt->IFG2_oneline[i + Nmax+1]);
         }
         
         P_IFG2[j] = fabs(creal(cpow(kdisc[j], 3) * f13_IFG2[j] * Pbin[j]));
     }
    
     
     gsl_interp_accel *acc_PIFG2
     = gsl_interp_accel_alloc();
     gsl_spline *spline_PIFG2
     = gsl_spline_alloc (gsl_interp_cspline, Nmax);
     gsl_spline_init(spline_PIFG2, kdisc, P_IFG2, Nmax);
     
     
     for (index_k=0; index_k < pnlpt->k_size; index_k++){
         
         pk_IFG2[index_k] = gsl_spline_eval(spline_PIFG2, pnlpt->k[index_k], acc_PIFG2);
     }
     
     gsl_spline_free (spline_PIFG2);
     gsl_interp_accel_free (acc_PIFG2);
     
     free(f13_IFG2);
     free(P_IFG2);
     
free(js);
free(kdisc);
free(Pbin);
free(Pdisc);
     
free(etam2);
free(cmsym2);
     
return _SUCCESS_;
}


// when you define the arguments of the function, don't forget to also change $nonlinear.h$ file

int nonlinear_pt_p13UV_integrate(
                                struct nonlinear_pt *pnlpt,
                                struct background *pba,
                                double * integrand_array,
                                int integrand_size,
                                int ia_size,
                                int index_ia_k,
                                int index_ia_pk,
                                int index_ia_sum,
                                int index_ia_ddsum,
                                double * sum
                                ) {
    
    double k,pk,integrand;
    int index_k;
    double bnorm = 1./(6. * pow(_PI_,2)) * pow(1 * pba->h,0);
    
    for (index_k=0; index_k < integrand_size; index_k++) {
        k = integrand_array[index_k*ia_size + index_ia_k];
        pk = integrand_array[index_k*ia_size + index_ia_pk];
        
        integrand = pk * bnorm;
      //  printf("%e  %e\n",k,pk);
        integrand_array[index_k*ia_size + index_ia_sum] = integrand;
    }
    
    /* fill in second derivatives */
    class_call(array_spline(integrand_array,
                            ia_size,
                            integrand_size,
                            index_ia_k,
                            index_ia_sum,
                            index_ia_ddsum,
                            _SPLINE_NATURAL_,
                            pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    /* integrate */
    class_call(array_integrate_all_spline(integrand_array,
                                          ia_size,
                                          integrand_size,
                                          index_ia_k,
                                          index_ia_sum,
                                          index_ia_ddsum,
                                          sum,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    return _SUCCESS_;
    
}

