#!/bin/sh 
    git add .
    git status
    git commit -m " code for 2D to 1D slab model "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
