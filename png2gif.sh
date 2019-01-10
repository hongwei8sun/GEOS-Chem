#!/bin/bash

rm -rf *.gif

convert  -delay  0  0_xy.png   -loop  0  0.gif

for ((i=1; i<=65; i++))

do

convert  -delay  0  ${i}_xy.png   -loop  0  ${i}.gif
convert  $[i-1].gif  ${i}.gif  ${i}.gif
rm -rf $[i-1].gif

done

mv  $[i-1].gif final.gif
