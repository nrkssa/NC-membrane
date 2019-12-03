#!/bin/bash
for i in CC_SOURCE CC_INCLUDE FC_SOURCE FC_INCLUDE
do
	cd $i
	sync-folder all `pwd`
	cd ../
done

