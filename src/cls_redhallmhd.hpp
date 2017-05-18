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
 * FILE: Definition of class "redhallmhd"
 *
 * DESCRIPTION: For defining and implementing the reduced mhd Hall physics for
 *              coronos. this class is responsible for "filling" lcsolve's data
 *              structures with the appropriate values - based on its
 *              "knowledge" of the physical model of * the plasma. These values
 *              are needed by lcsolve so that lcsolve can update its stack.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef CLS_REDHALLMHD
#define CLS_REDHALLMHD


#include <config.h>
#include "mpi.h"
#include "cls_parameter_map.hpp"
#include "cls_stack.hpp"
#include "cls_fft.hpp"
#include "nsp_constants.hpp"

#include <fstream>
#include <complex>
#include <vector>
#include <cstddef>
#include<iomanip>

#ifdef HAVE_CUDA_H
  #include "cls_redhallmhd_cuda_ext.hpp"
#else
#include <assert.h>
#endif

class redhallmhd
{

  private:

  ComplexArray roldlb;
  ComplexArray rnewlb;

  ComplexArray roldub;
  ComplexArray rnewub;

  void init_physics_data(               stack& run);

  void initTimeInc(                     stack& run);

/* ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ */

  void initU(                           stack& run);                  /* ~ U initialization functions                ~ */
  void computeFourierU(                 stack& run);

/* ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ */
  void computeRealU(                    stack& run);
  void initGauss(                       stack& run);
  void readUData(                       stack& run);

  void pLinzEnv(                        stack& run);
  void initBoundaries(                  stack& run);
  void countModes(                      stack& run);
  void initFootPointDriving(            stack& run);
  void initNoDrive(                     stack& run);

/* ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ */
  void retrieveOJ (                     stack& run);
/* ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ */
  void OfromP(                          stack& run );                 /* ~ Obtain vorticity from P                   ~ */
  void HfromA(                          stack& run );                 /* ~ Obtain H from A                           ~ */
//void JfromA(                          stack& run );                 /* ~ Obtain J from A                           ~ */


  void applyFootPointDrivingBC(std::string str_stp, stack& run );     /* ~ "pevol"                                   ~ */
  void applyLineTiedBC( std::string str_step, stack& run );           /* ~ pbot and p(:,n3) set to zero              ~ */

  void finalizeBoundaries(              stack& run );
  void finalizeFootPointDriving(        stack& run );                 /* ~                                           ~ */
  void finalizeLineTiedBoundaries(      stack& run );                 /* ~                                           ~ */

  public:

  parameter_map physics_data;

  fft fftw;

  ComplexArray P;                                                     /* ~ Stream function Phi in Fourier Space      ~ */
  ComplexArray O;                                                     /* ~  Vorticity storage  in Fourier Space      ~ */
  ComplexArray A;                                                     /* ~ flux function A in Fourier Space          ~ */
  ComplexArray J;                                                     /* ~ current density in Fourier Space          ~ */


  RealArray valfven;                                                  /* ~ Needed for Inhomogeneous RMHD             ~ */
  RealArray dvalfdz;
  RealArray nofz;
  RealArray dndz;
  RealArray umean;
  RealArray dudz;

  RealArray Elln;
  RealArray EllA;
  RealArray EllB;

  RealArray kpm;                                                      /* ~ kpm = k_p^-                               ~ */
  RealArray kpp;                                                      /* ~ kpp = k_p^+                               ~ */
  RealArray kmm;                                                      /* ~ kmm = k_m^-                               ~ */
  RealArray kmp;                                                      /* ~ kmp = k_m^+                               ~ */

  RealArray maxU;                                                     /* ~ for time-step determination               ~ */

  Real2DArray QtyVsZ;
  Real3DArray SpcVsZ;

  RealArray   ke;
  RealVar     dk;
  RealVar     dk_m1;
  RealVar     kb;
  RealVar     kf;

  int         isp;
  int         ikb;
  int         ikf;
  int         nk;

  void updatePAOJ( std::string str_step, stack& run );
/* ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ */
  void checkState( int pair, stack &run, std::string roc); 
/* ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ */
  void applyBC(    std::string str_step, stack& run );                 /* ~ Apply Boundary Conditions at current step ~ */
  void updateTimeInc(                    stack& run );

  void PfromO(                           stack& run );                 /* ~ Obtain P from vorticity                   ~ */
  void AfromH(                           stack& run );                 /* ~ Obtain A from H                           ~ */

  void evalElls(                         stack& run );                 /* ~ calculate l's and h's at each layer       ~ */
  void evalValf(                         stack& run );                 /* ~ calculate Va at each layer                ~ */
  void evalUmean(                        stack& run );                 /* ~ calculate Va at each layer                ~ */

                                                                      /* ~ Energy etc time-series related            ~ */

  void trackEnergies(    RealVar t_cur, stack& run );                 /* ~ update energy quantities between steps    ~ */
  void reportEnergyQs(                  stack& run );

  double evalTotalKineticEnergy(        stack& run );
  double evalTotalMagneticEnergy(       stack& run );
  double evalTotalVorticitySqd(         stack& run );
  double evalTotalCurrentSqd(           stack& run );
  double evalTotalGradCurrentSqd(       stack& run );
  double evalTotalFootPointKE(          stack& run );                 /* ~ Misnomer?                                 ~ */
  double evalTotalPoyntingFlux(         stack& run );                 /* ~ Poynting Flux                             ~ */
  double evalTotalHelicalEnergy(        stack& run );                 /* ~ volumetric energy contribution due to     ~ */
                                                                      /* ~ density gradient                          ~ */

                                                                      /* ~ Qty's vs z related                        ~ */
  void trackQtyVsZ(      RealVar t_cur, stack& run );                
  void reportQtyVsZ(                    stack& run );
                                                                      /* ~ Power Spectra related                     ~ */
  void trackPowerSpectra(RealVar t_cur, stack& run );
  void reportPowerSpectra(              stack& run );


  void physicsFinalize(                 stack& run );                 /* ~ end of subrun bookkeeping                 ~ */
                                                                      /* ~ this is just a stub right now, but will   ~ */
                                                                      /* ~ eventually make it possible to privatize  ~ */
                                                                      /* ~ some of the currently public procedures   ~ */

  redhallmhd();                                                       /* ~ Constructor (default)                     ~ */

/* ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ */
  redhallmhd(                           stack& run );                 /* ~ Constructor                               ~ */
/* ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ */

  ~redhallmhd();                                                      /* ~ Destructor                                ~ */


};

#endif
