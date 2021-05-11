#!/bin/bash

## Specifies the interpreting shell for the job.
#$ -S /bin/bash

## Specifies that all environment variables active within the qsub utility be exported to the context of the job.
#$ -V

## Specifies the parallel environment
#$ -pe smp 4

## Execute the job from the current working directory.
#$ -cwd 

## The  name  of  the  job.
#$ -N static30-2T

##send an email when the job ends
#$ -m e

##email addrees notification
#$ -M jl16@alumnes.udl.cat

##Passes an environment variable to the job
#$ -v  OMP_NUM_THREADS=2

## In this line you have to write the command that will execute your application.
./Static30-2T /share/apps/files/convolution/images/im01.ppm /share/apps/files/convolution/kernel/kernel3x3_Edge.txt /state/partition1/static30-2T.pm >result/static30-2T.txt