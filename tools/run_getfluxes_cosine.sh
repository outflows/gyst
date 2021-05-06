#!/bin/sh

# Give the job a name
#$ -N gyst_getfluxes

# set the shell
#$ -S /bin/sh

# set working directory on all host to
# directory where the job was started
#$ -cwd

# send all process STDOUT (fd 2) to this file
#$ -o job_output.txt

# send all process STDERR (fd 3) to this file
#$ -e job_output.err

# email information
#$ -m e

# Just change the email address. You will be emailed when the job has finished.
#$ -M rodrgust@oregonstate.edu

# generic parallel environment with X cores requested
#$ -pe orte 1

# Load a module, if needed
python/anaconda3-5.0.0.1

# Commands
python getenergy.py
