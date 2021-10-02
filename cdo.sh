for dir in /n/home12/hongwei/HONGWEI/data/GEOS-Chem/gcgrid/data/ExtData/GEOS_0.5x0.625/MERRA2/2020/*/  
do 
    dir=${dir%*/}
    echo ${dir##*/}
    echo 1
    pwd
    echo 3

    files=$(ls $dir/*.A3dyn.05x0625.nc4)
    echo ${files}

    cdo ensmean ${files} ${dir##*/}_ave.nc4

    popd
done
