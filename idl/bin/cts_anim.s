#!/bin/bash
 
################################################################################
#                                                                              #
# bash-shell script: cts_idl.s - Timothy J. Dennis 12 April 2013               #
#                             tdennis@gi.alaska.edu                            #
#                                                                              #
# Developed for use with the family of Heliospheric Physics codes originating  #
# with the work of C.S. Ng and including 'rmct_np' and 'fourfields.'           #
#                                                                              #
# Description: This script is designed to be executed from within a pbs script #
#              though it can be run from the command line if the user knows    #
#              which arguments to pass to it, and their proper ordering.       #
#              This is, however, not recommended for jobs requiring large      #
#              numbers of plots - particularly when the data is of high        #
#              resolution.                                                     #
#                                                                              #
#              It's function is to convert the eps files generated from the    #
#              coordinated execution of the c-shell script 'cts_idl.s', the    #
#              IDL batch script 'cts_batch', and the IDL procedure script -    #
#              and its dependencies - referenced therein to a set of gif       #
#              files. These gif files are then used to construct an animation  #
#              showing the time-evolution of a "horizontal" layer of the       #
#              solution space of the code from which the data was obtained.    #
#                                                                              #
#              As currently written the code can handle output for the fields  #
#              listed below:                                                   #
#                                                                              #
#              Field                                 cont_fld  Model Eq.'s     #
#                                                                              #
#              phi - 'Stream Function'               'p'       RMHD/FF         #
#              A   - 'Flux Function'                 'a'       RMHD/FF         #
#              J   - 'Current Density (z-comp)       'j'       RMHD/FF         #
#                                                                              #
#              Z   - 'B-field fluctations (z-comp)   'bz'      FF              #
#              V   - 'Velocity fluctuations (z-comp) 'vz'      FF              #
#                                                                              #
#                                                                              #
#              Also as currently written, this script should reside in a sub-  #
#              directory of the code run-directory called 'cts' though this    #
#              will likely change.                                             #
#                                                                              #
################################################################################

         qty=$1      # field/quantity to be contoured - see above
     res_str=$2      # string giving the x (or y) X z-resolution of the data
  desc_label=$3      # descriptive label for run
       n_slc=$4      # contoured layer

     cur_dir=`pwd`   # directory originating script invocation

cd cts

     cts_eps=`ls -1 $qty/eps/$qty'_ctr_'$res_str*.eps`

if [ ! -e $pfx/png ] # create gif directory if needed
  then 
    mkdir $qty/png
fi

       trunc=$qty'/eps/'

echo "trunc = " trunc

for plot in $cts_eps; do
        plot=${plot#$trunc}
    len_plot=${#plot}                                    # Get the file's basename and create a gif version
           n=`expr $len_plot - 4`                        # with convert
   base_name=${plot:0:$n}

   echo "creating " $base_name'.png'

   convert -density 300x300 $qty/eps/$base_name.eps -background  white -flatten -resize 1024X1024 $qty/png/$base_name.png
done

           cts_png=`ls $qty/png/$qty'_ctr_'$res_str*.png`
      cts_gif_anim=$qty'_cts_'$res_str'_'$desc_label'-anim.gif'
      cts_mpf_anim=$qty'_cts_'$res_str'_'$desc_label'-anim.mp4'

convert -delay 20 -dispose previous -loop 0 $cts_png $cts_gif_anim

$HOME/local/bin/ffmpeg -framerate 24 -pattern_type glob -i $qty/png/"*.png" -pix_fmt yuv420p $cts_mpf_anim

exit 0
