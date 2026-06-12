#!/bin/sh
#@$-q batch     # Submit to pipe queue "batch"
#@$-lM 15GB     # Request 1000MW memory
#@$-lT 15:30:00  # Request 2 1/2 hours CPU time
#@$-mb -me      # Send mail at beginning and end of request execution
#@$-eo          # Direct output from stderr to the stdout file
#@$-c 12
#@$             # REQUIRED! This line indicates the end of NQS options

# Shell script begins here.

cd $QSUB_WORKDIR         # cd back to the directory from
                         #   which the script was submitted
setenv F_PROGINF DETAIL  # request performance information
setenv F_TRACE YES
setenv F_FILEINF YES
