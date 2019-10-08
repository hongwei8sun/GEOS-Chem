#!/bin/sh 
    git add .
    git status
    git commit -m "Replacing GEOS-Chem input U/V data with modified new U/V data, so the U/V around polar point wouldn't be all the same value"
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
