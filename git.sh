#!/bin/sh 
    git add .
    git status
    git commit -m " rewrite plume model (polar coordinat) in Python "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
