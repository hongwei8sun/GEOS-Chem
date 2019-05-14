#!/bin/bash

rm -rf *.gif

convert  -delay  40  0_xy2.png   -loop  0  0_2.gif

for ((i=1; i<=30; i++))

do
echo $i
convert  -delay 40  ${i}_xy2.png   -loop  0  ${i}_2.gif
convert  -delay 40 $[i-1]_2.gif  ${i}_2.gif  ${i}_2.gif
rm -rf $[i-1]_2.gif

done

