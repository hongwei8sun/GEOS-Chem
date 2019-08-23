#!/bin/sh 
    git add .
    git status
    git commit -m " set z=0.0 in gaussian analytical results !!! "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
