#!/bin/bash
cd $1
for i in $(ls -d AB*)
do
	echo "renaming  $i "
	dir1=`echo $i |cut -f 1 -d '('`
	dir2=`echo $i |cut -f 1 -d ")"|cut -f 2 -d "(" `
	echo $dir1 $dir2
	mkdir -pv $dir1
	mv -v $i $dir1/$dir2
	cd $dir1/$dir2
	for j in {1..3}
	do
		k=$((j+3))
		mv -v ENS-$j ../../../Equm/zstart-0/$dir1/$dir2/ENS-$k
	done
	cd ../../
done

