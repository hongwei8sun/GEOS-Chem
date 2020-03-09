#!/bin/sh 
    git add .
    git status
    git commit -m " dissolve the plume when there is only 1 ring left "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
