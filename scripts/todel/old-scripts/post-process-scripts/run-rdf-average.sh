#!/bin/bash
script_dir=`pwd`
script_name='Make-rdf-Histogram.py'
cd $1
for d1 in blen*/
do
	cd $d1
	for d2 in Flex*/
	do
		cd $d2
		for d3 in */
		do
			cd $d3
			for d4 in zstart*/
			do
				cd $d4
				for d5 in kappa*/
				do
					cd $d5
					for d6 in AB*/
					do
						cd $d6
						cp -vf $script_dir/$script_name .
						python $script_name
						rm -vf $script_name
						cd ../
					done
					cd ../ #d5
				done
				cd ../ #d4
			done
			cd ../ #d3
		done
		cd .. #d2 
	done
	cd ..  #d1
done

