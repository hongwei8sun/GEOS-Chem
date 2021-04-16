#!/bin/sh 
    git add .
    git status
    git commit -m "1 plume in 60 degree"
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
