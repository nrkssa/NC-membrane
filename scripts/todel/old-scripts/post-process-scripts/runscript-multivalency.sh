#!/bin/bash
script_name='Make-Histogram.py'
script_dir='/media/5f0becfe-1890-4202-afc5-e6325625e84c/Nanocarriers-code/Fluctuating-Surface/Data-Fluct-Surface/'
source_file='Multivalency-Data.dat'
cd $1
for i in $(ls -d zstart*)
do
	cd $i
	for j in $(ls -d AB*)
	do
		cd $j
		for k in $(ls -d N-*)
		do
			cd $k
			cp -vf $script_dir/$script_name .
			python $script_name $source_file
			rm -vf $script_name
			cd ..
		done
		cd ..
	done
	cd ..
done

