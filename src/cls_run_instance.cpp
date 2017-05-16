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
 *        FILE: Implementation of class "run_instance"
 *
 * DESCRIPTION: To be added prior to public release
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "cls_run_instance.hpp"

/* ~ Constructors ~ */

run_instance::run_instance() {

  int is_init;
  int mpi_stat;

  mpi_stat = MPI_Initialized(&is_init);

  if (is_init) {

    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    std::string pname;
    std::string padjust;

    int tot_proc;

    std::string run_start_time = getTime();
    std::string node_name      = getNode();
    mpi_stat                   = MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);

    int ndevice; int * ptr_ndevice; ptr_ndevice = &ndevice;

#ifdef HAVE_CUDA_H
    ri_cuext.getDeviceCount(ptr_ndevice);
#else
    ndevice = 0;
#endif

    padjust.assign("rfx");

    pname.assign("run_start_time");
    run_data.emplace(pname, run_start_time, padjust);

    pname.assign("node_name");
    run_data.emplace(pname, node_name,      padjust);

    pname.assign("rank");
    run_data.emplace(pname, rank,           padjust);

    pname.assign("tot_proc");
    run_data.emplace(pname, tot_proc,       padjust);

    pname.assign("ndevice");
    run_data.emplace(pname, ndevice,        padjust);

  }
  else {

      mpi_stat = MPI_Abort(MPI_COMM_WORLD, MPI_ERR_REQUEST);

  }
}

/* ~ Destructor ~ */

run_instance::~run_instance() {

}

/* ~ utilities ~ */

std::string run_instance::getTime() {

  time_t utc_sec    = time(0);
  tm    *gm_time    = gmtime(&utc_sec);
  char  *c_str_time = asctime(gm_time);

  std::string str_time(c_str_time);

  return str_time;

}

std::string run_instance::getNode() {

  const int          maxlen = 132;
        char         nodename[maxlen];
        int          success;
        std:: string str_node;

  success = gethostname(nodename, maxlen);

  if (success == 0) {

    str_node.assign(nodename);

  }
  else {

    str_node.assign("run_instance: WARNING - getNode failed");

  }
  return str_node;
}
