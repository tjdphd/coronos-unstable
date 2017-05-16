#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#

  AC_PREREQ([2.63])
  AC_INIT([coronos], [0.01], [tdennis10@alaska.edu])
  AC_LANG(C++)
  AM_INIT_AUTOMAKE([foreign])
  AC_CONFIG_SRCDIR([src/coronos.cpp])
  AC_CONFIG_HEADERS([config.h])
  AC_CONFIG_MACRO_DIR([m4])

# Get the architecture

# -> if ARCH is 64 are lib64's gauranteed to exist? 

  AC_MSG_CHECKING([architecture])
            ARCH=`uname -m`
  if test x$ARCH==xx86_64
  then
          SUFFIX="64"
  else
          SUFFIX=""
  fi
  AC_MSG_RESULT([$ARCH])
  AC_MSG_NOTICE([setting suffix to $SUFFIX])

# Establish precision to be used 

  AH_TEMPLATE([LD_PRECISION_H], [Define to one if using long double precision])
  AH_TEMPLATE([OD_PRECISION_H], [Define to one if using ordinary double precision])

  AC_ARG_WITH([long_doubles], AS_HELP_STRING([--with-long_doubles],
                                     [use long double precision (default is no)])],
                                     [],[])

  if test -n "$with_long_doubles"
  then
    AC_DEFINE([LD_PRECISION_H])
  else
    AC_DEFINE([OD_PRECISION_H])
  fi

# get C++ compiler/mpi wrapper
    
  CRS_CXX="mpic++ mpiCC CC "
  AC_PROG_CXX([$CRS_CXX])

  AX_COMPILER_VERSION

  AC_MSG_NOTICE([found compiler vendor $ax_cv_cxx_compiler_vendor])
  AC_MSG_NOTICE([found compiler version $ax_cv_cxx_compiler_version])

  PORTLAND=`expr match "$ax_cv_cxx_compiler_vendor" '\(portland\)'`
     INTEL=`expr match "$ax_cv_cxx_compiler_vendor" '\(intel\)'`
       GNU=`expr match "$ax_cv_cxx_compiler_vendor" '\(gnu\)'`

  if  test -n "$PORTLAND"
  then VENDOR=$PORTLAND
  fi
  if  test -n "$INTEL"
  then VENDOR=$INTEL
  fi
  if  test -n "$GNU"
  then VENDOR=$GNU
  fi

# Find mpi.h

  AC_MSG_CHECKING([where to find MPI_PATH])

  if ! test -n "$MPI_PATH"
  then
    if  test -n "$MPIHOME"
    then
        MPI_PATH="$MPIHOME"
    else
      if test -n "$MPI_HOME"
      then
        MPI_PATH=$MPI_HOME
      else
        AC_MSG_ERROR([cannot find path to Message Passing Interface library. Consider setting MPI_PATH manually],[1])
      fi
    fi
  fi

  MPI_PATH_INC="$MPI_PATH/include"
  if test -n "$MPI_PATH"
  then
    AC_CHECK_HEADERS([$MPI_PATH_INC/"mpi.h"],[],[],[])
  else
   AC_MSG_ERROR([cannot locate mpi.h],[1])
  fi
  
  MPI_PATH_LIB="$MPI_PATH/lib$SUFFIX"
  if  [ ! test -e "$MPI_PATH_LIB"  ||  ! test -d "$MPI_PATH_LIB" ]
  then
    MPI_PATH_LIB="$MPI_PATH/lib"
  fi

  AC_MSG_NOTICE([MPI_PATH_LIB = $MPI_PATH_LIB])

# locate FFTW path

  FFTW_H="fftw3.h"
  FFTW_LIBS="-lm -lfftw3"

  AC_MSG_CHECKING([whether FFTW_PATH is set])

  if ! test -n "$FFTW_PATH"
  then
    AC_MSG_RESULT([no, searching...])
    if test -n "$FFTWHOME"
    then
    FFTW_PATH=$FFTWHOME
    else 
      if test -n "$FFTW_ROOT"

      then
      FFTW_PATH=$FFTW_ROOT
      else
        AC_MSG_ERROR([cannot find path to FFTW library. Check if fftw module is loaded, or consider setting FFTW_PATH manually],[1])
      fi
    fi
  else
    AC_MSG_RESULT([yes])
  fi

  AC_MSG_NOTICE([FFTW_PATH = $FFTW_PATH])

  FFTW_PATH_INC=$FFTW_PATH/include
  AC_CHECK_HEADERS([$FFTW_PATH_INC/$FFTW_H],[
  AC_MSG_NOTICE([please disregard the superfluous warning above regarding rejection by the pre-processor.])
  ],[],[])

  FFTW_PATH_LIB="$FFTW_PATH/lib$SUFFIX"
  if  [ ! test -e "$FFTW_PATH_LIB" || ! test -d "$FFTW_PATH_LIB" ]
  then
    FFTW_PATH_LIB="$FFTW_PATH/lib"
  fi

  AC_MSG_NOTICE([FFTW_PATH_LIB = $FFTW_PATH_LIB])

# Prepare to test for CUDA

  AH_TEMPLATE([DONT_HAVE_CUDA_H], [Define to one if not using cuda])
  AH_TEMPLATE(     [HAVE_CUDA_H], [Define to one if using cuda])

  AC_ARG_WITH([cuda], AS_HELP_STRING([--with-cuda=PATH],
                                     [use gpu acceleration (default is no)])],
                                     [],[])    


  if test -n "$with_cuda"
  then 
    AC_MSG_CHECKING([whether nvcc is available])
    AC_CHECK_PROG([NVCC],[nvcc],[nvcc],[""],[],[])
    if test -n "$NVCC"
    then
      AC_MSG_RESULT([yes])
      AC_MSG_CHECKING([seeking path to CUDA resources])

      if ! test -n "$CUDA_PATH"
      then
        if  test -n "$CUDAHOME"
        then
          CUDA_PATH="$CUDAHOME"
        else
          AC_DEFINE([DONT_HAVE_CUDA_H])
          AC_MSG_ERROR([Cannot find path to CUDA resources, consider setting CUDA_PATH manually.][1])
        fi
      fi

      if test -n "$CUDA_PATH"
      then
            CUDA_LIBS=" -lcudart"
               CUDA_H="cuda.h"
        CUDA_PATH_INC="$CUDA_PATH/include"
        CUDA_PATH_LIB="$CUDA_PATH/lib$SUFFIX"
        if  [ ! test -e "$CUDA_PATH_LIB" || ! test -d "$CUDA_PATH_LIB" ]
        then
          CUDA_PATH_LIB="$FFTW_PATH/lib"
        fi
        AC_CHECK_HEADERS([$CUDA_PATH_INC/$CUDA_H],[],[],[])

        # still need cufft and cublas

           CUFFT_H="cufft.h"
          CUFFTW_H="cufftw.h"
        CUFFT_LIBS=" -lcufft -lcufftw"

        AC_MSG_CHECKING([header $CUFFT_H])
        AC_CHECK_HEADERS([$CUDA_PATH_INC/$CUFFT_H],[],[],[])
        AC_MSG_CHECKING([header $CUFFTW_H])
        AC_CHECK_HEADERS([$CUDA_PATH_INC/$CUFFTW_H],[],[],[])

        CUBLAS_LIBS=" -lcublas"
        CUBLAS_H="cublas_v2.h"
        AC_MSG_CHECKING([header $CUBLAS_H])
        AC_CHECK_HEADERS([$CUDA_PATH_INC/$CUBLAS_H],[],[],[])

        # if we get this far we have everything needed for GPU acceleration

        AC_DEFINE([HAVE_CUDA_H])
CORONOS_SYSTEM='cls_run_instance_cuda_ext.cu.$(OBJEXT) cls_canvas_cuda_ext.cu.$(OBJEXT) cls_stack_cuda_ext.cu.$(OBJEXT) cls_lcsolve_cuda_ext.cu.$(OBJEXT) cls_redhallmhd_cuda_ext.cu.$(OBJEXT) cls_fft_cuda_ext.cu.$(OBJEXT)'
        echo "CORONOS_SYSTEM = " $CORONOS_SYSTEM
        AC_SUBST([CORONOS_SYSTEM])
        AC_SUBST([NVCC])

      else
        AC_DEFINE([DONT_HAVE_CUDA_H])
        AC_MSG_ERROR([Path to CUDA resources is not set, consider setting CUDA_PATH manually.][1])
      fi
    else

      AC_DEFINE([DONT_HAVE_CUDA_H])
      AC_MSG_ERROR([nvcc CUDA compiler wrapper not detected. Load cuda module and try again.],[1])

    fi
  else
    AC_DEFINE([DONT_HAVE_CUDA_H])
    # maybe still need to do non-cuda stuff here
  fi
  
# Wrap up

   CRS_CPPFLAGS="-I$MPI_PATH_INC -I$CUDA_PATH_INC -I$FFTW_PATH_INC -I.. "
    CRS_LDFLAGS="-L$CUDA_PATH_LIB -L$FFTW_PATH_LIB -L$MPI_PATH_LIB "
       CRS_LIBS="$CUDA_LIBS $CUFFT_LIBS  $CUBLAS_LIBS  $FFTW_LIBS"

  
  AC_SUBST(CRS_CPPFLAGS)
  AC_SUBST(CRS_LIBS)
  AC_SUBST(CRS_LDFLAGS)

  AC_CONFIG_FILES([Makefile
                 src/Makefile])
  AC_OUTPUT
