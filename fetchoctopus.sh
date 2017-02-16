#!/bin/bash
# Simple script that fetches data from scratch to the current working directory.
# Expects the presence of file job.log created by r.abin.

if [[ ! -e *.log ]];then
   echo "ERROR: no file .log does not exist. Exiting now..."
   exit 1
fi

NODE=$(head -3 job.log)
JOB=$(tail -4 job.log)
KAM=$(tail -5 job.log)

# copy all data from scratch if it is newer (-u switch)
# and preserve the timestamps (-p switch)
ssh -n $NODE "cp -r -u -p $JOB/* $KAM"

