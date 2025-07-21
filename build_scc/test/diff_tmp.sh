#!/bin/bash

for f in openacc/tmp/* 
do
  #f = ref
  #g = file name
  g=$(basename $f)
  echo "==> $f <=="
  diff $f tmp/$g
done
