#!/bin/sh 
    git add .
    git status
    git commit -m " add runge kutta into lagrangian (assume wind speed doesn't change with time, but change with location) "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
