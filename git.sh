#!/bin/sh 
    git add .
    git status
    git commit -m " Conventional GEOS-Chem: aircraft releases aerosols one by one "
#    git commit -m " aircraft releases aerosols one by one "
#    git commit -m " Lagrangian figure: aircraft releases aerosols "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
