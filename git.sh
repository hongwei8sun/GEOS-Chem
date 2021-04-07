#!/bin/sh 
    git add .
    git status
    git commit -m "Aircraft_plume_evolution.gif"
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
