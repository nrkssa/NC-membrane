#!/bin/bash
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
			for l in $(ls -d ENS-*)
			do
				echo $l
				rm -vf $l/SYS_STATE/MEMB/*.vtu
				rm -vf $l/SYS_STATE/NC-DATA/*.vtu
				rm -vf $l/SYS_STATE/DUMP/*.dump
			done
			cd ..
		done
		cd ..
	done
	cd ..
done

