#!/bin/bash

mkdir copy_out
mkdir copy_out/summary_pngs
mkdir copy_out/bead_pngs
mkdir copy_out/full_bead_segment_pngs
mkdir copy_out/quant
cp *_log.txt copy_out/
cp *.csv copy_out/quant/
cp *_D_*.png copy_out/bead_pngs/
cp *_B_*.png copy_out/summary_pngs/
cp *_C_*.png copy_out/full_bead_segment_pngs

