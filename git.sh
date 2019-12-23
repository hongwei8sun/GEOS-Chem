#!/bin/sh 
    git add .
    git status
    git commit -m " fix the error about particles concentrated in South Polar "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
