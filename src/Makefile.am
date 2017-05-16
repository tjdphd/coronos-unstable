%.cu.o:	%.cu
	$(NVCC) -arch=sm_20 $(CRS_CPPFLAGS) -o $@ -c $<
bin_PROGRAMS     = coronos
coronos_SOURCES  = coronos.cpp cls_run_instance.cpp  \
			             cls_canvas.cpp  \
				      cls_stack.cpp  \
				    cls_lcsolve.cpp  \
			         cls_redhallmhd.cpp  \
			      cls_parameter_map.cpp  \
			          cls_parameter.cpp  \
                                        cls_fft.cpp
                                                       

EXTRA_coronos_SOURCES = cls_run_instance_cuda_ext.cu \
                        cls_canvas_cuda_ext.cu       \
                        cls_stack_cuda_ext.cu        \
                        cls_lcsolve_cuda_ext.cu      \
                        cls_redhallmhd_cuda_ext.cu   \
                        cls_fft_cuda_ext.cu

coronos_CPPFLAGS      = $(CRS_CPPFLAGS)
coronos_LDFLAGS       = $(CRS_LDFLAGS)

coronos_LDADD         = $(CRS_LIBS) $(CRS_LDFLAGS) $(CORONOS_SYSTEM)
coronos_DEPENDENCIES  = $(CORONOS_SYSTEM)

