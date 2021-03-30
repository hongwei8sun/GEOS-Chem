#!/bin/sh 
    git add .
    git status
    git commit -m " instead of injecting plume in the whole simulation, plume now can be injected only in a period at the beginning"
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
