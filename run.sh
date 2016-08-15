#!/bin/bash

python run.py four data/ --out_pdf output_figs/all.pdf --output_prefix output --L 600  --sw_window_size 2000

#rm figs/ambig.pdf
#python run.py four /data/mtsinai/2016_05_13_GR37_HG002_hapcalls/ambig_calls/ figs/ambig --L 600;

#rm figs/hetero.pdf;
#python run.py four /data/mtsinai/2016_05_13_GR37_HG002_hapcalls/heterozygous_hap1_calls/ figs/hetero.pdf --L 600;

#rm figs/homo.py;
#python run.py four /data/mtsinai/2016_05_13_GR37_HG002_hapcalls/homozygous_calls/ figs/hetero.pdf --L 600;

#python run.py sixteen data/all_files/ figs/all.pdf
#python run.py classify data/ambig_calls/ figs/all.pdf
