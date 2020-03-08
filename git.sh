#!/bin/sh 
    git add .
    git status
    git commit -m " add: minus the extra amount of injected aerosol when dissolve the plume "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
