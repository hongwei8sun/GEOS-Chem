#!/bin/sh 
    git add .
    git status
    git commit -m "aircraft release aerosol boxes between (-30N,30N) for 1 year"
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
