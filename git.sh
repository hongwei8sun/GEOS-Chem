#!/bin/sh 
    git add .
    git status
    git commit -m " compare combining rings used in model with analytical results "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
