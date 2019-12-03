#!/bin/bash
filename1=equm-runscript-1.py
filename2=equm-runscript-2.py
n=`ps -ef | grep equm- | wc -l`
while [ $n -ge 5 ]; do sleep 300;  n=`ps -ef | grep equm- | wc -l`; done
python $filename1
send-email $filename1 ramn@seas.upenn.edu
sleep 20

n=`ps -ef | grep equm- | wc -l`
while [ $n -ge 5 ]; do sleep 300;  n=`ps -ef | grep equm- | wc -l`; done
python $filename2
send-email $filename2 ramn@seas.upenn.edu
