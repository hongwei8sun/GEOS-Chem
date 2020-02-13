#!/bin/sh 
    git add .
    git status
    git commit -m " the version before optimizing to speed up"
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
