#include "cls_run_instance_cuda_ext.hpp"

run_instance_cuda_ext::run_instance_cuda_ext() { }

void run_instance_cuda_ext::getDeviceCount(int * devices ) { cudaError_t err = cudaGetDeviceCount(devices); }

run_instance_cuda_ext::~run_instance_cuda_ext() { }
