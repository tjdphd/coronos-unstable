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
 *        FILE: CUDA extension for definition of class "fft"
 *
 * DESCRIPTION: To be added prior to public release
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <config.h>
#include <cuda.h>
#include <cufft.h>
#include <cufftw.h>
#include "cls_stack.hpp"

class fft_cuda_ext
{

  private:
  
  public:

    cufftHandle     cu_p_lay_for;                             /* ~ For establishing plans for forward   ~ */
    cufftHandle     cu_p_lay_rev;                             /* ~ reverse FFT's of layers              ~ */

    cufftDoubleComplex *cu_cplx_out;
    cufftDoubleReal    *cu_r_in;

    ComplexVar * host_cplx_out;
    RealVar    * host_real_in;

    void cufftwInitialize( stack& run);
    void cufftwFinalize();

    void cufftwForwardIC(RealArray&    Rin, ComplexArray& Cout);
    void cufftwReverseIC(ComplexArray& Cin, RealArray&    Rout);

    fft_cuda_ext();
   ~fft_cuda_ext();

};
