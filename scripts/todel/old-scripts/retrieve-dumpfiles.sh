#!/bin/bash
curdir=`pwd`
cd $1
dir1=`echo $1 | cut -f 2-3 -d '_'|cut -f 1 -d '/'`
echo $dir1

	for i1 in $(ls -d blen-*)
	do
	sf=`echo $i1 | cut -f 4 -d '-'`
	fdir1=$curdir/memb-conf/$dir1/'sf-'$sf
	n=`ls $fdir/*.dump | wc -l`
	echo $fdir $n
	cd $i1
		for i2 in $(ls -d Flex*)
		do
		cd $i2
			for i3 in $(ls -d *)
			do
			cd $i3
				for i4 in $(ls -d zstart*)
				do
				cd $i4
					for i5 in $(ls -d kappa*)
					do
					cd $i5
						fdir=$fdir1/$i5
						for i6 in $(ls -d AB*)
						do
						cd $i6
							for i7 in $(ls -d ENS*)
							do
							cd $i7
								cd SYS_STATE/DUMP
								n1=`ls *.dump | wc -l`
								n1=$((n1-1))
								if find $fdir/ -maxdepth 0 -empty |read
								then
								n=1
								else
								n=`ls $fdir/*.dump | wc -l`
								n=$((n+1))
								fi
								cp -rv ./membrane_state-$n1.dump $fdir/membrane_state-$n.dump
								cd ../../
	
							cd ../
							done
						cd ../
						done
					cd ../
					done
				cd ../
				done
				
			cd ..
			done
		cd ../
		done
	cd ../
	done
cd ../
