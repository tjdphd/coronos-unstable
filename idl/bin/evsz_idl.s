#!/bin/csh

#################################################################################
#                                                                               #
# c-shell script: evsz_idl.s - Timothy J. Dennis 10 April 2014                  #
#                             tdennis@gi.alaska.edu                             #
#                                                                               #
# Developed for use with the family of Heliospheric Physics codes originating   #
# with the work of C.S. Ng and including 'rmct_np' and 'fourfields.'            #
#                                                                               #
# Description: This script is designed to be executed from within a pbs script  #
#              though it can be run from the command line if the user knows     #
#              which arguments to pass to it, and their proper ordering.        #
#              This is, however, not recommended for jobs requiring large       #
#              numbers of plots - particularly when the data is of high         #
#              resolution.                                                      #
#                                                                               #
#              It's function is to prepare a compute node's environment for     #
#              submission of  the batch script 'evsz_batch' to an instantiation #
#              of idl on a compute node and allows for the generation of large  #
#              numbers of plots based on high-resolution data to be             #
#              produced 'automatedly' so as to avoid the need to stress the     #
#              resources of a login node, or the patience of the code-jockey    #
#              (that would me or you) who needs said contour plots to move his  #
#              research forward.                                                #
#                                                                               #
# Usage Example:                                                                #
#                                                                               #
#              ./evsz_idl.s ts 0 99                                             #
#                                                                               #
#              Note: the argument values given here must be consistent with     #
#                    similar arguments chosen in the scripts used when          #
#                    initiating the run that produced the data to be post-      #
#                    processed.                                                 #
#                                                                               #
#################################################################################
#
# Arguments to 'cts_idl.s' are ultimately to be passed to 'cts_batch'. This is
# accomplished by defining corresponding environment variables which may
# subsequently be accessed by IDL via the function 'GETENV'
#
############################
                           #
setenv EVSZ_DSC_LAB     $1 # - string coding for information about the data
setenv EVSZ_RES_STR     $2 #
setenv EVSZ_EFLD        $3 # - string coding for information about the data
setenv EVSZ_FIRST_STEP  $4 # - First data-set in series to be plotted
setenv EVSZ_LAST_STEP   $5 # - Last data-set in series to be plotted
setenv EVSZ_MINMAX_MODE $6 # - either 'global' or 'local'
                           #
############################
#
#
# Next we need to be sure that the appropriate sub-directory
# that is assumed to exist by the batch script will in fact
# exist.
#
############################################
                                           #
if( ! -d evsz) then                        # upon first run the output
  mkdir evsz                               # directory will not yet exist
  mkdir -p evsz/$EVSZ_EFLD/eps             # directory will not yet exist
else if( ! -d evsz/$EVSZ_EFLD) then        #
  mkdir -p evsz/$EVSZ_EFLD/eps             # directory will not yet exist
else if( ! -d evsz/$EVSZ_EFLD/eps) then    #
  mkdir evsz/$EVSZ_EFLD/eps                # directory will not yet exist
endif                                      #
                                           #
############################################
#
# Next we load the idl module, tell it where to find IDL sources and
# tell it the name of the batch script we want it to execute
#
#########################################################################################
                                                                                        #
source /etc/profile.d/modules.csh                                                       #
module load idl                                                                         # load idl onto the compute node
                                                                                        #
setenv EVSZ_SRCDIR $PWD/idl/pro/evsz                                                    #
setenv IDL_PATH    $EVSZ_SRCDIR                                                         #  
setenv IDL_PATH    $IDL_DIR/lib:$IDL_PATH                                               #            
setenv IDL_PATH    $IDL_DIR/lib/graphics:$IDL_PATH                                      #                     
setenv IDL_PATH    $IDL_DIR/lib/bridges:$IDL_PATH                                       #                    
setenv IDL_PATH    $IDL_DIR/lib/datatypes:$IDL_PATH                                     #                      
setenv IDL_PATH    $IDL_DIR/lib/dicomex:$IDL_PATH                                       #                    
setenv IDL_PATH    $IDL_DIR/lib/enterprise:$IDL_PATH                                    #                       
setenv IDL_PATH    $IDL_DIR/lib/enterprise/ese:$IDL_PATH                                #                           
setenv IDL_PATH    $IDL_DIR/lib/hook:$IDL_PATH                                          #                 
setenv IDL_PATH    $IDL_DIR/lib/imsl:$IDL_PATH                                          #                 
setenv IDL_PATH    $IDL_DIR/lib/obsolete:$IDL_PATH                                      #                     
setenv IDL_PATH    $IDL_DIR/lib/utilities:$IDL_PATH                                     #                      
setenv IDL_PATH    $IDL_DIR/lib/wavelet:$IDL_PATH                                       #                    
setenv IDL_PATH    $IDL_DIR/lib/wavelet/bitmaps:$IDL_PATH                               #                            
setenv IDL_PATH    $IDL_DIR/lib/wavelet/source:$IDL_PATH                                #                           
setenv IDL_PATH    $IDL_DIR/lib/wavelet/data:$IDL_PATH                                  #                         
setenv IDL_PATH    $IDL_DIR/lib/itools/:$IDL_PATH                                       #
setenv IDL_PATH    $IDL_DIR/lib/itools/framework:$IDL_PATH                              #                             
setenv IDL_PATH    $IDL_DIR/lib/itools/components:$IDL_PATH                             #                              
setenv IDL_PATH    $IDL_DIR/lib/itools/ui_widgets:$IDL_PATH                             #                              
setenv IDL_PATH    $IDL_DIR/bin:$IDL_PATH                                               #            
setenv IDL_PATH    $IDL_DIR/bin.linux.x86_64:$IDL_PATH                                  #                         
setenv IDL_PATH    $IDL_DIR/bin.linux.x86_64/dm:$IDL_PATH                               #                            
setenv IDL_PATH    $IDL_DIR/bin.linux.x86_64/dm/bin:$IDL_PATH                           #                                
setenv IDL_PATH    $IDL_DIR/bin.linux.x86_64/dm/demo:$IDL_PATH                          #                                 
setenv IDL_PATH    $IDL_DIR/bin.linux.x86_64/dm/help:$IDL_PATH                          #                                 
setenv IDL_PATH    $IDL_DIR/bin.linux.x86_64/dm/lib:$IDL_PATH                           #                                
setenv IDL_PATH    $IDL_DIR/bin.linux.x86_64/dm/locale:$IDL_PATH                        #                                   
setenv IDL_PATH    $IDL_DIR/bin.linux.x86_64/dm/locale/en_US:$IDL_PATH                  #
setenv IDL_PATH    $IDL_DIR/bin.linux.x86_64/dm/locale/en_US/LC_MESSAGES:$IDL_PATH      #
setenv IDL_PATH    $IDL_DIR/bin/xml:$IDL_PATH                                           #
setenv IDL_PATH    $IDL_DIR/bin/make_rt:$IDL_PATH                                       #
                                                                                        #
setenv IDL_STARTUP $PWD/idl/bin/evsz_batch                                              # idl batch script to execute
                                                                                        #
idl                                                                                     # invoke idl
                                                                                        #
#########################################################################################
#
# Now clean things up.
#
####################################
                                   #
module unload idl                  #
unsetenv IDL_STARTUP               #
unsetenv IDL_PATH                  #
unsetenv EVSZ_DSC_LAB              #
unsetenv EVSZ_RES_STR              #
unsetenv EVSZ_EFLD                 #
unsetenv EVSZ_FIRST_STEP           #
unsetenv EVSZ_LAST_STEP            #
unsetenv EVSZ_MINMAX_MODE          #
unsetenv EVSZ_SRCDIR               #
unsetenv IDL_PATH                  #
unsetenv IDL_STARTUP               #
                                   #
exit 0                             #
                                   #
####################################
##
###
#### The End ####
