#!/bin/sh 
    git add .
    git status
    git commit -m " Analytical skewed plume concentration distribution "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
