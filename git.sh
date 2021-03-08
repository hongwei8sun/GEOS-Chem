#!/bin/sh 
    git add .
    git status
    git commit -m " set plume length to 2000m, change Kchem value "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
