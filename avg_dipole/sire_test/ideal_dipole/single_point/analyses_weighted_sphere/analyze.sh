#!/bin/bash


for f in l_* ; do 
    echo $f
    bondlength=${f##*_}
    python average_distance.py "$f/*_deg/computed.dat" "$f/*_deg/theory.dat" $bondlength
    wait
done
