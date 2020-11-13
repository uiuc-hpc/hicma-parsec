#!/bin/bash

for i in `ls | grep run`; do
  echo $i
  sbatch $i
  sleep 1
done
