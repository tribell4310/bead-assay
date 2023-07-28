#!/bin/bash

for filename in *; do
  echo $filename
  if [ "${filename: -5}" == ".tiff" ]
  then
    prefix=$(echo $filename | cut -d'.' -f 1)
    suffix="_%d.tiff"
    magick convert $filename "$prefix$suffix"
    rm $filename
  fi
done

