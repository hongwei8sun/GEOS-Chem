#!/bin/sh 
    git add .
    git status
    git commit -m " update lagrange_mod.F90"
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
