#!/bin/bash

#@ job_type = serial

#@ class = Small

#@ environment = COPY_ALL

#@ queue

#@ output = job.out

#@ error = job.err

 

wmake > C.out





