#!/bin/sh 
    git add .
    git status
    git commit -m "Sigma 50m, 100m, 150m"
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
