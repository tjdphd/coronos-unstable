#!/bin/bash

 desc_label=$1
    res_str=$2
        qty=$3
 first_step=$4
  last_step=$5

zero_str='000000000000000000000000000000'
wild_str='??????????????????????????????'
     sep='_'
     hph='-'
 cur_dir=`pwd` # directory originating script invocation

cd evsz

len_last_step=`expr length $last_step`
#
if [ $len_last_step -lt 3 ]
then
       len_dif=`expr 3 - $len_last_step`
           pos=`expr 30 - $len_dif`
 trnc_zero_str=${zero_str:$pos}
     last_step=$trnc_zero_str$last_step
fi


if [ $len_last_step -lt 3 ]
then
  trnc_wild_str='???'
else
            pos=`expr 30 - $len_last_step`
  trnc_wild_str=${wild_str:$pos}
fi
      evsz_eps=`ls -1 $qty/eps/'evsz'$sep$qty$sep$res_str$sep'zl'$sep*$hph$trnc_wild_str'.eps'`
      evsz_gif_anim=$qty'_evsz'$hph'steps'$sep$first_step$hph$last_step$sep$res_str$hph$desc_label$hph'anim.gif'
      evsz_mpf_anim=$qty'_evsz'$hph'steps'$sep$first_step$hph$last_step$sep$res_str$hph$desc_label$hph'anim.mp4'

if  ! [   -e $qty/png ]
then
  mkdir $qty/png
fi


trunc=$qty'/eps/'
echo "trunc = " $trunc
for plot in $evsz_eps; do
        plot=${plot#$trunc}
    len_plot=${#plot}                                    # Get the file's basename and create a gif version
           n=`expr $len_plot - 4`                        # with convert
   base_name=${plot:0:$n}
   echo "creating = " $base_name'.png'
   convert -density 300x300 $qty/eps/$base_name.eps -background white -flatten -resize 1024x1024 $qty/png/$base_name.png
done
#
evsz_png=`ls -1 $qty/png/'evsz'$sep$qty$sep$res_str$sep'zl'$sep*$hph$trnc_wild_str'.png'`

convert -delay 20 -dispose previous -loop 0 $evsz_png $evsz_gif_anim

$HOME/local/bin/ffmpeg -framerate 24 -pattern_type glob -i $qty/png/"*.png" -pix_fmt yuv420p $evsz_mpf_anim

exit 0
