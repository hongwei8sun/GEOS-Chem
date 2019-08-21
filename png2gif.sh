#!/bin/bash

rm -rf *.gif

convert  -delay  60  0_xy.png   -loop  0  0.gif

for ((i=1; i<=24; i++))

do
echo $i
convert  -delay  60  ${i}_xy.png   -loop  0  ${i}.gif
convert -delay 40  $[i-1].gif  ${i}.gif  ${i}.gif
rm -rf $[i-1].gif

done

