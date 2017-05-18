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
 *        FILE: Definition of class "lcsolve"
 *
 * DESCRIPTION: solver class for longcope solver designed to work on canvas 
 *              class type "stack"
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef CLS_LCSOLVE
#define CLS_LCSOLVE

#include <config.h>
#include "mpi.h"
#include<iomanip>
#include <fstream>
#include <complex>
#include <vector>
#include <cstddef>
#include "nsp_constants.hpp"
#include "cls_stack.hpp"

using namespace constants;
#include "cls_redhallmhd.hpp"

#ifdef HAVE_CUDA_H
  #include "cls_lcsolve_cuda_ext.hpp"
#else
  #include <assert.h>
#endif


class lcsolve
{

   private:

   RealArray    SE0; /* ~ same for different layers of U's      ~ */
   RealArray    SE1;
   RealArray    SE2;
   RealArray    SE3;

   RealArray    SI0; /* ~ same for different layers of U's      ~ */
   RealArray    SI1;
   RealArray    SI2;
   RealArray    SI3;

   ComplexArray B0; /* ~ different for different layers of U's ~ */
   ComplexArray B1;
   ComplexArray B2;
   ComplexArray B3;

   ComplexArray D0; /* ~ different for different layers of U's ~ */
   ComplexArray D1;
   ComplexArray D2;
   ComplexArray D3;

   ComplexArray A0; /* ~ different for different layers of U's ~ */
   ComplexArray A1;
   ComplexArray A2;
   ComplexArray A3;

   void    createFields(stack& run );
   void    destroyFields();

   RealVar maxdU(RealArray& dx, RealArray&  dy, int i_grid, int i_layers);

   void    averageAcrossLayers(        stack& run, int shift_sign, RealArray&    dx, RealArray&  dy);
   void    averageAcrossLayers(        stack& run, int shift_sign, ComplexArray& dx                );

   void    Step( std::string str_step, stack& run );

   void    setS(    std::string str_step,   stack& run, redhallmhd& physics );
   void    setB(    std::string str_step,   stack& run, redhallmhd& physics );
   void    setD(    std::string str_step,   stack& run, redhallmhd& physics );
   void    setAi(                           stack& run, redhallmhd& physics );

/* ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ */

   void    partialsInXandY(                 stack& run, redhallmhd& physics, ComplexArray& U, RealArray& Ux, RealArray& Uy);

   void    bracket(                         stack& run, redhallmhd& physics, 
                                                                             ComplexArray& BrKt, 
                                                                                RealArray& dx1, 
                                                                                RealArray& dy1, 
                                                                                RealArray& dx2, 
                                                                                RealArray& dy2
                  );

/* ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ */

   public:

// lcsolve();                                    /* ~ Constructors                         ~ */

    lcsolve(   stack& run );
   ~lcsolve();                                   /* ~ Destructor                           ~ */

#ifdef HAVE_CUDA_H
  lcsolve_cuda_ext lcs_cuext;
#endif

/* ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ */

   void Loop(                                     stack& run ); /* ~ stepping and such                     ~ */

/* ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ REQUIRES FFT ~ */

   void passAdjacentLayers( std::string str_step, stack& run );

};

#endif
