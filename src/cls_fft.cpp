/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 *
 * CORONOS||SONOROC - Version 0.1
 *
 * (S)ynthesized  (O)bject-based (N)umerical (O)bservatory for (R)HMHD [also RMHD and IRHMD] with (O)ptional (C)UDA-acceleration
 *
 * AUTHOR: Timothy J. Dennis
 *         tdennis10@alaska.edu
 *
 * CONTRIBUTORS:
 *
 *         C. S. Ng
 *         LiWei Lin
 *         Others to be included prior to public release
 *
 * copyright 2014-2016 
 *
 * Space Physics and Aeronomy
 * Geophysical Institute
 * University of Alaska, Fairbanks
 *
 * All Rights Reserved.
 *
 * This version of the code is pre-public release.
 * Please contact the author if you are not certain
 * you have an up-to-date working copy.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*
 *        FILE: Implementation of class "fft"
 *
 * DESCRIPTION: To be added prior to public release
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "cls_fft.hpp"

#ifdef DONT_HAVE_CUDA_H

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~ non-CUDA Fourier-related ~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwInitialize( stack& run ) {

  fftwKInit(  run );
  fftwrtInit( run );


  int n1;    run.stack_data.fetch("n1",   &n1   );
  int n2;    run.stack_data.fetch("n2",   &n2   );
  int nr_in; run.stack_data.fetch("n1n2", &nr_in);
  int nc_out;

  nc_out    = n1 * (((int)(0.5*n2)) + 1);

/* ~ Forward field/layer transforms: ~ */

#ifdef LD_PRECISION_H
  r_in      = (RealVar *)               fftwl_malloc(sizeof(RealVar)               * nr_in  );
  cplx_out  = (ComplexVar *) fftwl_malloc(sizeof(ComplexVar) * nc_out );
  p_lay_for = fftwl_plan_dft_r2c_2d(n1, n2, r_in, reinterpret_cast<fftwl_complex*>(cplx_out), FFTW_EXHAUSTIVE);
#elif defined OD_PRECISION_H
  r_in      = (RealVar *)               fftw_malloc(sizeof(RealVar)               * nr_in  );
  cplx_out  = (ComplexVar *) fftw_malloc(sizeof(ComplexVar) * nc_out );
  p_lay_for = fftw_plan_dft_r2c_2d(n1, n2, r_in, reinterpret_cast<fftw_complex*>(cplx_out), FFTW_EXHAUSTIVE);
#endif

/* ~ Reverse field/layer transforms: ~ */

#ifdef LD_PRECISION_H
  cplx_in   = (ComplexVar *) fftwl_malloc(sizeof(ComplexVar) * nc_out );
  r_out     = (RealVar *)    fftwl_malloc(sizeof(RealVar)    * nr_in  );
  p_lay_rev = fftwl_plan_dft_c2r_2d(n1, n2, reinterpret_cast<fftwl_complex*>(cplx_in), r_out, FFTW_EXHAUSTIVE);
#elif defined OD_PRECISION_H
  cplx_in   = (ComplexVar *) fftw_malloc(sizeof(ComplexVar) * nc_out );
  r_out     = (RealVar *)    fftw_malloc(sizeof(RealVar)    * nr_in  );
  p_lay_rev = fftw_plan_dft_c2r_2d(n1, n2, reinterpret_cast<fftw_complex*>(cplx_in), r_out, FFTW_EXHAUSTIVE);
#endif

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwFinalize() {

#ifdef LD_PRECISION_H
 fftwl_destroy_plan(p_lay_for);
 fftwl_destroy_plan(p_lay_rev);

 fftwl_free(r_in);
 fftwl_free(r_out);

 fftwl_free(cplx_in); 
 fftwl_free(cplx_out); 
#elif defined OD_PRECISION_H
 fftw_destroy_plan(p_lay_for);
 fftw_destroy_plan(p_lay_rev);

 fftw_free(r_in);
 fftw_free(r_out);

 fftw_free(cplx_in); 
 fftw_free(cplx_out); 
#endif

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwKInit(stack& run) {               /* ~ initialize the wave number arrays kx, ky, k2, and inv_k2 ~ */

  RealArray& kx     = run.kx; 
  RealArray& ky     = run.ky; 
  RealArray& k2     = run.k2; 
  RealArray& inv_k2 = run.inv_k2;

  int n1; run.stack_data.fetch("n1",    &n1);    /* ~ number of coordinates in x                               ~ */
  int n2; run.stack_data.fetch("n2",    &n2);    /* ~ number of coordinates in y                               ~ */ 
  int nc; run.stack_data.fetch("n1n2c", &nc);    /* ~ number of Fourier space points per layer                 ~ */

  int n1h    = (((int)(half*n1)) + 1);
  int n2h    = (((int)(half*n2)) + 1);

  RealArray kp;                                  /* ~ temporary definition of convenience                      ~ */

  kp.assign(n1,zero);

      kx.assign(nc,zero);                        /* ~ these are all lcstack members and carry the wave-number  ~ */
      ky.assign(nc,zero);                        /* ~ information                                              ~ */
      k2.assign(nc,zero);
  inv_k2.assign(nc,zero);

                                                 /* ~ the initialization as done in gpu port of old code       ~ */
  int ndx;
  int i_cpy;

  for (int i = 0;   i < n1h; i++ ) kp[i] = ((RealVar) i) * two_pi;
  for (int i = n1h; i < n1;  i++ ) kp[i] = -kp[n1 - i];

  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2h; j++) {
   
      ndx     = (i * n2h) + j;
      kx[ndx] = kp[j];
      ky[ndx] = kp[i];

    }
  }

//for (int i = 0; i < n1; i++) {
//  for (int j = 0; j < n1h; j++) {

//    ndx     = (i * n2h) + j;
//    k2[ndx] = (kx[ndx] * kx[ndx]) +  (ky[ndx] * ky[ndx]);

//    if (std::abs(k2[ndx]) > teensy) inv_k2[ndx] = one / k2[ndx];
//    else inv_k2[ndx] = zero;

//  }
//}

//for (int l = 0; l < nc; ++l) {

//      if (l  < (nc - (n2/2)-1)) {
//        if ((l == 0) || (l %((n2/2)+1) !=0) ){
//          ; // do nothing
//        }
//        else {
//               i_cpy    = (l/((n2/2)+1));
//               if ( i_cpy < n2/2 ) {
//                 kx[l] = kx[i_cpy];
//                 ky[l] = ky[i_cpy];
//               }
//        }
//      }
//      else {
//               i_cpy    = (l - nc + n2/2+2);
//             if (i_cpy < 0 || i_cpy > n2/2+1) {
//               std::cout << "fftwKInit: WARNING - i_cpy = "             << i_cpy << " for l = " << l <<  std::endl;
//             }
//             else {

//                 kx[l] = kx[i_cpy];
//                 ky[l] = ky[i_cpy];
//             }
//      }
//}

//kp.resize(0);

/* ~ TEST ~  TEST ~  TEST ~  TEST ~  TEST ~  TEST ~  TEST ~  TEST ~  TEST ~  TEST ~  TEST ~  TEST ~ */

  for ( int l = 0; l < nc; ++ l) { k2[l] = (kx[l] * kx[l]) + (ky[l] * ky[l]); 
                                   if (k2[l] > teensy) { inv_k2[l] = one /  k2[l]; }
                                   else                { inv_k2[l] = zero;         }
                                 }
//for (int i = 0; i < n1; i++) {
//  for (int j = 0; j < n1h; j++) {

//    ndx     = (i * n2h) + j;
//    k2[ndx] = (kx[ndx] * kx[ndx]) + (ky[ndx] * ky[ndx]);


//    if (std::abs(k2[ndx]) > teensy) inv_k2[ndx] = one / k2[ndx];
//    else inv_k2[ndx] = zero;
//    else inv_k2[ndx] = huge;


//  }
//}
/* ~ TEST ~  TEST ~  TEST ~  TEST ~  TEST ~  TEST ~  TEST ~  TEST ~  TEST ~  TEST ~  TEST ~  TEST ~ */
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwKFree(stack& run ) {                            /* ~ also possibly a bit excessive ~ */

    RealArray&     kx = run.kx;
    RealArray&     ky = run.ky;
    RealArray&     k2 = run.k2;
    RealArray& inv_k2 = run.inv_k2;

        kx.resize(0);
        ky.resize(0);
        k2.resize(0);

    inv_k2.resize(0);
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwrtInit(stack& run) {             /* ~ initialize de-aliasing array              ~ */

  int rank;     MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  RealArray& kx = run.kx;
  RealArray& ky = run.ky;

  int n1;    run.stack_data.fetch("n1",    &n1);    /* ~ number of x-coordinates                   ~ */
  int n2;    run.stack_data.fetch("n2",    &n2);    /* ~ number of y-coordinates                   ~ */
  int n1n2c; run.stack_data.fetch("n1n2c", &n1n2c); /* ~ number of Fourier space points per layer  ~ */

  RealArray::size_type nc;                          /* ~ a vector size-type version of n1n2c       ~ */
  nc  = n1n2c;

  int n1h = (((int)(half*n1)) + 1);
  int n2h = (((int)(half*n2)) + 1);

  RealVar kx_max = two_pi * half * n1;
  RealVar ky_max = two_pi * half * n2;
  RealVar k_max;

  if (kx_max >= ky_max) { k_max = ky_max;}
  else                  { k_max = kx_max;}

  RealVar threshold = two_thirds * k_max;

  rt.assign(nc, zero);                              /* ~ create space in rt                        ~ */

  RealArray::size_type ndx;

  for (unsigned i = 0; i < n1; i++) {               /* ~ initialization                            ~ */
    for (unsigned j = 0; j < n1h; j++) {

      ndx = i * n2h + j;

      if (  ky[ndx] == zero) { 
        if ( std::abs(kx[ndx]) < k_max) {
          rt[ndx] = one;
        }
        else {
          rt[ndx] = zero;
//        rt[ndx] = one;
        }
      }
      else if ( std::abs(ky[ndx]) >  threshold  || std::abs(kx[ndx]) >  threshold  ) {
        rt[ndx] = zero;
//      rt[ndx] = one;
      }
      else {
        rt[ndx] = one;
      }

  }
 }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwrtFree() {

 rt.resize(0);

}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwForwardAll( stack& run ) {


  InputOutputArray& U   = run.U;               /* ~ raw input array                         ~ */

  ComplexArray& U0      = run.U0;
  ComplexArray& U1      = run.U1;
  ComplexArray& U2      = run.U2;
  ComplexArray& U3      = run.U3;

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */
  int srun;     run.palette.fetch("srun", &srun);
  int rank;     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */
  int n1;       run.stack_data.fetch("n1"   , &n1      );
  int n2;       run.stack_data.fetch("n2"   , &n2      );
  int n1n2c;    run.stack_data.fetch("n1n2c", &n1n2c   );
  int n1n2;     run.stack_data.fetch("n1n2" , &n1n2    );
  int n_layers; run.stack_data.fetch("iu2"  , &n_layers);
  int n_flds;   run.stack_data.fetch("iu3"  , &n_flds  );

  RealVar scale    = (RealVar) one / ((RealVar) (n1n2));

  ComplexArray::size_type nc = n1n2c;

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */

//if ((srun > 1) && (rank == 0)) {
//  for (unsigned k = 0; k < n1n2; k++) {
//    if (U[k][0][0] != zero) { std::cout << "U[" << k << "][0][0] = " << U[k][0][0] << std::endl; }
//  }
//}

/* ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ TEST ~ */

  unsigned strt_idx;

  for (int i_f = 0; i_f < n_flds; i_f++)  {
    for ( int i_l = 0; i_l < n_layers; i_l++) {
      for (unsigned k = 0; k < n1n2; k++) { r_in[k]     = zero; }
      for (unsigned k = 0; k < n1n2; k++) { r_in[k]     = U[k][i_l][i_f]; }

#ifdef LD_PRECISION_H
      fftwl_execute(p_lay_for);
#elif defined OD_PRECISION_H
      fftw_execute(p_lay_for);
#endif

      strt_idx = i_l * nc;

      switch(i_f) {
      case(0) :
        for (unsigned k = 0; k < nc;   k++) { U0[(strt_idx + k)] = scale*cplx_out[k]*rt[k]; }
        U0[(strt_idx)].real()=zero;
        break;
      case(1) :  
        for (unsigned k = 0; k < nc;   k++) { U1[(strt_idx + k)] = scale*cplx_out[k]*rt[k]; }
        U1[(strt_idx)].real()=zero;
        break;
      case(2) :  
        for (unsigned k = 0; k < nc;   k++) { U2[(strt_idx + k)] = scale*cplx_out[k]*rt[k]; }
        U2[(strt_idx)].real()=zero;
        break;
      case(3) :  
        for (unsigned k = 0; k < nc;   k++) { U3[(strt_idx + k)] = scale*cplx_out[k]*rt[k]; }
        U3[(strt_idx)].real()=zero;
        break;
      }
    }
  }


}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwReverseAll( stack& run ) {


  InputOutputArray& U  = run.U;               /* ~ raw input array                         ~ */

  ComplexArray&     U0 = run.U0;
  ComplexArray&     U1 = run.U1;
  ComplexArray&     U2 = run.U2;
  ComplexArray&     U3 = run.U3;

  int n1;       run.stack_data.fetch("n1"   , &n1      );
  int n2;       run.stack_data.fetch("n2"   , &n2      );
  int n1n2c;    run.stack_data.fetch("n1n2c", &n1n2c   );
  int n1n2;     run.stack_data.fetch("n1n2" , &n1n2    );
  int n_layers; run.stack_data.fetch("iu2"  , &n_layers);
  int n_flds;   run.stack_data.fetch("iu3"  , &n_flds  );

  ComplexArray::size_type nc = n1n2c;

  unsigned strt_idx;

  for (int i_f = 0; i_f < n_flds; i_f++)  {
    for ( int i_l = 0; i_l < n_layers; i_l++) {
      
      strt_idx = (i_l * nc);

      for (unsigned   k = 0; k < nc; k++) { cplx_in[k] = (ComplexVar) zero; }
      switch(i_f) {
      case(0) :
        U0[strt_idx].real()=zero;
        for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U0[strt_idx + k];  }
        break;
      case(1) :
        U1[strt_idx].real()=zero;
        for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U1[strt_idx + k];  }
        break;
      case(2) :
        U2[strt_idx].real()=zero;
        for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U2[strt_idx + k];  }
        break;
      case(3) :
        U3[strt_idx].real()=zero;
        for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U3[strt_idx + k];  }
        break;
      }

      for (int k = 0; k < n1n2; k++) { r_out[k] = zero; }
#ifdef LD_PRECISION_H
      fftwl_execute(p_lay_rev);
#elif defined OD_PRECISION_H
      fftw_execute(p_lay_rev);
#endif

      for (int k = 0; k < n1n2; k++) { U[k][i_l][i_f] =  r_out[k]; }

    }
  }


}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwForwardAll( stack& run, ComplexArray& Oout, ComplexArray& Jout ) {


  InputOutputArray& AUX  = run.AUX;             /* ~ raw input array                         ~ */

  int n1;       run.stack_data.fetch("n1"   , &n1      );
  int n2;       run.stack_data.fetch("n2"   , &n2      );
  int n1n2c;    run.stack_data.fetch("n1n2c", &n1n2c   );
  int n1n2;     run.stack_data.fetch("n1n2" , &n1n2    );
  int n_layers; run.stack_data.fetch("iu2"  , &n_layers);

  RealVar scale    = (RealVar) one/((RealVar) (n1n2));

  ComplexArray::size_type nc = n1n2c;

  unsigned strt_idx;

  for (int i_f = 0; i_f < 2; i_f++)  {
    for ( int i_l = 0; i_l < n_layers; i_l++) {

      for (unsigned k = 0; k < n1n2; k++) { r_in[k]     = zero; }
      for (unsigned k = 0; k < n1n2; k++) { r_in[k]     = AUX[k][i_l][i_f]; }

#ifdef LD_PRECISION_H
      fftwl_execute(p_lay_for);
#elif defined OD_PRECISION_H
      fftw_execute(p_lay_for);
#endif

      strt_idx = i_l * nc;

      switch(i_f) {
      case(0) :
        for (unsigned k = 0; k < nc;   k++) { Oout[(strt_idx + k)] = scale*cplx_out[k]*rt[k]; }
        Oout[strt_idx].real()=0;
        break;
      case(1) :  
        Jout[strt_idx].real()=0;
        for (unsigned k = 0; k < nc;   k++) { Jout[(strt_idx + k)] = scale*cplx_out[k]*rt[k]; }
        break;
      }
    }
  }


}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwReverseAll( stack& run, ComplexArray& Oin, ComplexArray& Jin ) {


  InputOutputArray& AUX = run.AUX;               /* ~ raw input array                         ~ */

  int n1;       run.stack_data.fetch("n1"   , &n1      );
  int n2;       run.stack_data.fetch("n2"   , &n2      );
  int n1n2c;    run.stack_data.fetch("n1n2c", &n1n2c   );
  int n1n2;     run.stack_data.fetch("n1n2" , &n1n2    );
  int n_layers; run.stack_data.fetch("iu2"  , &n_layers);

  ComplexArray::size_type nc = n1n2c;

  unsigned strt_idx;

  for (int i_f = 0; i_f < 2; i_f++)  {
    for ( int i_l = 0; i_l < n_layers; i_l++) {
      
      strt_idx = (i_l * nc);

      for (unsigned   k = 0; k < nc; k++) { cplx_in[k] = (ComplexVar) zero; }
      switch(i_f) {
      case(0) :
        Oin[strt_idx].real()=0;
        for (unsigned k = 0; k < nc; k++) { cplx_in[k] = Oin[strt_idx + k];  }
        break;
      case(1) :
        Jin[strt_idx].real()=0;
        for (unsigned k = 0; k < nc; k++) { cplx_in[k] = Jin[strt_idx + k];  }
        break;
      }

      for (int k = 0; k < n1n2; k++) { r_out[k] = zero; }
#ifdef LD_PRECISION_H
      fftwl_execute(p_lay_rev);
#elif defined OD_PRECISION_H
      fftw_execute(p_lay_rev);
#endif

      for (int k = 0; k < n1n2; k++) { AUX[k][i_l][i_f] =  r_out[k]; }

    }
  }


}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwForwardLayerofField ( stack& run, int i_l, int i_f ) {


  InputOutputArray& U = run.U;               /* ~ raw input array                         ~ */

  ComplexArray& U0 = run.U0;
  ComplexArray& U1 = run.U1;
  ComplexArray& U2 = run.U2;
  ComplexArray& U3 = run.U3;

  int n1; 
  run.stack_data.fetch("n1"   , &n1      );
  int n2; 
  run.stack_data.fetch("n2"   , &n2      );
  int n1n2c; 
  run.stack_data.fetch("n1n2c", &n1n2c   );
  int n1n2; 
  run.stack_data.fetch("n1n2" , &n1n2    );
  int n_layers; 
  run.stack_data.fetch("iu2"  , &n_layers);
  int n_flds;
  run.stack_data.fetch("iu3"  , &n_flds  );

  RealVar scale    = (RealVar) one/((RealVar) (n1n2));

  ComplexArray::size_type nc = n1n2c;

  for (unsigned k = 0; k < n1n2; k++) { r_in[k]     = zero; }
  for (unsigned k = 0; k < n1n2; k++) { r_in[k]     = U[k][i_l][i_f]; }
  for (unsigned k = 0; k < nc;   k++) { cplx_out[k] = (ComplexVar) zero; }

#ifdef LD_PRECISION_H
  fftwl_execute(p_lay_for);
#elif defined OD_PRECISION_H
  fftw_execute(p_lay_for);
#endif

  unsigned strt_idx = i_l * nc;

  switch(i_f) {
  case(0) :
    for (unsigned k = 0; k < nc;   k++) { U0[(strt_idx + k)] = scale*cplx_out[k]*rt[k]; }
    U0[strt_idx].real()=zero;
    break;
  case(1) :  
    U1[strt_idx].real()=zero;
    for (unsigned k = 0; k < nc;   k++) { U1[(strt_idx + k)] = scale*cplx_out[k]*rt[k]; }
    break;
  case(2) :  
    U2[strt_idx].real()=zero;
    for (unsigned k = 0; k < nc;   k++) { U2[(strt_idx + k)] = scale*cplx_out[k]*rt[k]; }
    break;
  case(3) :  
    U3[strt_idx].real()=zero;
    for (unsigned k = 0; k < nc;   k++) { U3[(strt_idx + k)] = scale*cplx_out[k]*rt[k]; }
    break;
  }


}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwReverseLayerofField ( stack& run, int i_l, int i_f) {


  InputOutputArray& U = run.U;               /* ~ raw input array                         ~ */

  ComplexArray&    U0 = run.U0;
  ComplexArray&    U1 = run.U1;
  ComplexArray&    U2 = run.U2;
  ComplexArray&    U3 = run.U3;

  int n1;       run.stack_data.fetch("n1"   , &n1      );
  int n2;       run.stack_data.fetch("n2"   , &n2      );
  int n1n2c;    run.stack_data.fetch("n1n2c", &n1n2c   );
  int n1n2;     run.stack_data.fetch("n1n2" , &n1n2    );
  int n_layers; run.stack_data.fetch("iu2"  , &n_layers);
  int n_flds;   run.stack_data.fetch("iu3"  , &n_flds  );

  ComplexArray::size_type nc = n1n2c;

  unsigned strt_idx = (i_l * nc);

  for (unsigned   k = 0; k < nc; k++) { cplx_in[k] = czero; }
  switch(i_f) {
  case(0) :
    U0[strt_idx].real()=zero;
    for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U0[strt_idx + k];  }
    break;
  case(1) :
    U1[strt_idx].real()=zero;
    for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U1[strt_idx + k];  }
    break;
  case(2) :
    U2[strt_idx].real()=zero;
    for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U2[strt_idx + k];  }
    break;
  case(3) :
    U3[strt_idx].real()=zero;
    for (unsigned k = 0; k < nc; k++) { cplx_in[k] = U3[strt_idx + k];  }
    break;
  }

  for (int k = 0; k < n1n2; k++) { r_out[k] = zero; }

#ifdef LD_PRECISION_H
  fftwl_execute(p_lay_rev);
#elif defined OD_PRECISION_H
  fftw_execute(p_lay_rev);
#endif

  for (int k = 0; k < n1n2; k++) { U[k][i_l][i_f] = r_out[k]; }


}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwForwardRaw( stack& run, RealArray& Rin, ComplexArray& Cout) {


  int n1n2c; run.stack_data.fetch("n1n2c", &n1n2c   );
  int n1n2;  run.stack_data.fetch("n1n2" , &n1n2    );
  int iu2;   run.stack_data.fetch("iu2"  , &iu2     );

  RealVar scale       = ((RealVar) one)/((RealVar) (n1n2));

  unsigned c_strt_idx = 0;
  unsigned r_strt_idx = 0;

  for (unsigned i_l   = 0; i_l < iu2; ++i_l) {

    c_strt_idx        = i_l * n1n2c;
    r_strt_idx        = i_l * n1n2;

    for (unsigned k   = 0 ; k < n1n2 ; ++k) { r_in[k]              =  zero;               }
    for (unsigned k   = 0 ; k < n1n2c; ++k) { cplx_out[k]          = czero;               }
    for (unsigned k   = 0 ; k < n1n2 ; ++k) { r_in[k]              = Rin[r_strt_idx + k]; }

#ifdef LD_PRECISION_H
    fftwl_execute(p_lay_for);
#elif defined OD_PRECISION_H
     fftw_execute(p_lay_for);
#endif

    for (unsigned k   = 0 ; k < n1n2c ; ++k){ Cout[c_strt_idx + k] = scale*cplx_out[k]*rt[k];}

    Cout[c_strt_idx].real() = zero;

  }


}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwReverseRaw( stack& run, ComplexArray& Cin, RealArray& Rout) {


  int n1n2c; run.stack_data.fetch("n1n2c", &n1n2c );
  int n1n2;  run.stack_data.fetch("n1n2" , &n1n2  );
  int iu2;   run.stack_data.fetch("iu2"  , &iu2   );

  unsigned c_strt_idx = 0;
  unsigned r_strt_idx = 0;

  for (unsigned i_l = 0; i_l < iu2; ++i_l) {

    c_strt_idx        = i_l * n1n2c;
    r_strt_idx        = i_l * n1n2;

    Cin[c_strt_idx].real() = zero;

    for (unsigned k = 0 ; k < n1n2c; ++k  ) { cplx_in[k]           = czero;               }
    for (unsigned k = 0 ; k < n1n2 ; ++k  ) { r_out[k]             =  zero;               }
    for (unsigned k = 0 ; k < n1n2c; ++k  ) { cplx_in[k]           = Cin[c_strt_idx + k]; }

#ifdef LD_PRECISION_H
    fftwl_execute(p_lay_rev);
#elif defined OD_PRECISION_H
    fftw_execute(p_lay_rev);
#endif

    for (unsigned k = 0; k < n1n2; ++k) { Rout[r_strt_idx + k] = r_out[k];}
  
  }


}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwReverseIC(ComplexArray& Cin, RealArray& Rout ) {


  int n1n2c; n1n2c   = Cin.size();
  int n1n2;  n1n2    = Rout.size();

  Cin[0].real()      = zero;

  for (unsigned k = 0 ; k < n1n2c; k++) { cplx_in[k]  = czero;        }
  for (unsigned k = 0 ; k < n1n2 ; k++) { r_out[k]    = zero;         }
  for (unsigned k = 0 ; k < n1n2c; k++) { cplx_in[k]  = Cin[k];       }

#ifdef LD_PRECISION_H
  fftwl_execute(p_lay_rev);
#elif defined OD_PRECISION_H
  fftw_execute(p_lay_rev);
#endif

  for (unsigned k = 0 ; k < n1n2 ; k++) {Rout[k]     = (r_out[k]);}


}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void fft::fftwForwardIC( RealArray& Rin, ComplexArray& Cout) {


  int n1n2c; n1n2c = Cout.size();
  int n1n2;  n1n2  = Rin.size();

  RealVar scale    = ((RealVar) one)/((RealVar) (n1n2));

  for (unsigned k = 0 ; k < n1n2 ; ++k) {r_in[k]     =  zero;       }
  for (unsigned k = 0 ; k < n1n2c; ++k) {cplx_out[k] = czero;       }
  for (unsigned k = 0 ; k < n1n2 ; ++k) {r_in[k]     = Rin[k];      }

#ifdef LD_PRECISION_H
  fftwl_execute(p_lay_for);
#elif defined OD_PRECISION_H
  fftw_execute(p_lay_for);
#endif

  for (unsigned k = 0 ; k < n1n2c; k++) {Cout[k]     = scale * cplx_out[k]*rt[k]; }

  Cout[0].real()   = zero;


}

#endif

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
