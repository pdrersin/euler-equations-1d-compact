#!/bin/bash

dx1=0.001
dx2=`echo $dx1*2|bc`
xdisc=0.7

x0=0
x1L=`echo $xdisc-$dx1|bc`
x1R=$xdisc
x2=1

echo "# cell-centered x-values" > points.dat
echo "`seq $x0 $dx1 $x1L`" >> points.dat
echo "`echo $xdisc-$dx1|bc`" >> points.dat
echo "`seq $x1R $dx2 $x2`" >> points.dat
