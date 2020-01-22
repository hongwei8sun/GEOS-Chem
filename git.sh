#!/bin/sh 
    git add .
    git status
    git commit -m " update flexchem.F90 "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
