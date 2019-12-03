#!/bin/bash
for i in $(seq 1 $1)
do
	dirname=SAMP-RUN$i
	rm -rfv $dirname
	mkdir -pv $dirname
	mkdir -pv $dirname/RUNDIR
	cp -rvf ./NC_MEMB $dirname/RUNDIR
	cp -rvf CC_INPUT_FILES FC_INPUT_FILES SYS_STATE $dirname
	cd $dirname/RUNDIR
	echo run >script
	gnome-terminal -e "gdb -x script NC_MEMB" &
	cd ../../
done
	
