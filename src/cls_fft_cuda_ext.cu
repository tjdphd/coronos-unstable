#include "cls_fft_cuda_ext.hpp"

fft_cuda_ext::fft_cuda_ext() {

}

void fft_cuda_ext::cufftwInitialize( stack& run) {

  cufftResult res_cufft;

  int n1;    run.stack_data.fetch("n1",   &n1   );
  int n2;    run.stack_data.fetch("n2",   &n2   );
  int n3;    run.stack_data.fetch("n3",   &n3   );

  int nr_in; run.stack_data.fetch("n1n2", &nr_in);
  int n1n2 ; run.stack_data.fetch("n1n2", &n1n2 );
  int n1n2c; run.stack_data.fetch("n1n2c", &n1n2c);
  int nc_out;

  nc_out    = n1 * (((int)(0.5*n2)) + 1);

  std::cout << "cufftwInitialize: nc_out = " << nc_out << std::endl;
  std::cout << "cufftwInitialize: n1n2c  = " << n1n2c  << std::endl;

  cudaMalloc((void**)&cu_cplx_out,sizeof(cufftDoubleComplex)*n1n2c);
  if ( cudaGetLastError() != cudaSuccess){ std::cout << "cufftwInitialize: unable to allocate cplx_out"     << std::endl;}

//  cudaHostAlloc( (void**) &host_cplx_out, sizeof(ComplexVar)*n1n2c,cudaHostAllocPortable );
//  cudaHostAlloc( (void**) &host_real_in,  sizeof(RealVar)   *n1n2c,cudaHostAllocPortable );

  host_cplx_out = (ComplexVar *)malloc(sizeof(ComplexVar)*n1n2c);
  host_real_in  = (RealVar    *)malloc(sizeof(RealVar)   *n1n2c);

  cudaMalloc((void**)&cu_r_in,    sizeof(cufftDoubleReal)   *n1n2 );
  if ( cudaGetLastError() != cudaSuccess){ std::cout << "cufftwInitialize: unable to allocate r_in"         << std::endl;}

  res_cufft = cufftPlan2d(&cu_p_lay_for, n1, n2, CUFFT_D2Z);
  if (res_cufft != CUFFT_SUCCESS) { std::cout        << "cufftwInitialize: could not create plan cu_p_lay_for" << std::endl;}
  res_cufft = cufftPlan2d(&cu_p_lay_rev, n1, n2, CUFFT_Z2D);
  if (res_cufft != CUFFT_SUCCESS) { std::cout        << "cufftwInitialize: could not create plan cu_p_lay_rev" << std::endl;}

}

void fft_cuda_ext::cufftwFinalize() {

  cudaError_t res_err;
  cufftResult res_cufft;

  res_cufft      = cufftDestroy(cu_p_lay_for);
  if (res_cufft != CUFFT_SUCCESS) {std::cout << "cufftwFinalize: could not destroy plan cu_p_lay_for" << std::endl;}
  res_cufft      = cufftDestroy(cu_p_lay_rev);
  if (res_cufft != CUFFT_SUCCESS) {std::cout << "cufftwFinalize: could not destroy plan cu_p_lay_rev" << std::endl;}

  res_err        = cudaFree(cu_cplx_out);
  if (res_err   != cudaSuccess) { std::cout  << "cufftwFinalize: unable to deallocate cplx_out"       << std::endl;}
  res_err        = cudaFree(cu_r_in);
  if (res_err   != cudaSuccess) { std::cout  << "cufftwFinalize: unable to deallocate r_in"           << std::endl;}

  free(host_cplx_out);
  free(host_real_in );

}

void fft_cuda_ext::cufftwForwardIC(RealArray& Rin, ComplexArray& Cout ) {

  int n1n2c; n1n2c = Cout.size();
  int n1n2;  n1n2  = Rin.size();

  RealVar scale    = ((RealVar) one)/((RealVar) (n1n2));

  for (unsigned k = 0 ; k < n1n2 ;  ++k) {host_real_in[k] =  Rin[k];                       }

  cudaMemcpy(                 cu_r_in, host_real_in, ( sizeof(RealVar)*n1n2),cudaMemcpyHostToDevice    );
  cufftExecD2Z( cu_p_lay_for, cu_r_in, cu_cplx_out);
  cudaMemcpy(host_cplx_out,            cu_cplx_out,  ( sizeof(ComplexVar)*n1n2c),cudaMemcpyDeviceToHost);

  for (unsigned k = 0 ; k < n1n2c ; ++k) {Cout[k]         = scale * host_cplx_out[k]; }

  /* ~ should be able to access rt to allow for dealiasing (need to friend?) ~ */ 

}

void fft_cuda_ext::cufftwReverseIC(ComplexArray& Cin, RealArray& Rin )  {

}

fft_cuda_ext::~fft_cuda_ext() {

}
