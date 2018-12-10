# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

alias ls='ls --color=tty'

export PATH=$HOME/mpich/bin:$PATH

# User specific aliases and functions

module load gcc/7.1.0-fasrc01
module load hdf5/1.10.1-fasrc03
module load zlib/1.2.8-fasrc09
module load szip/2.1-fasrc02
module load libpng/1.6.25-fasrc01
module load jpeg/6b-fasrc02
module load gcc/7.1.0-fasrc01 jasper/1.900.1-fasrc02
module load perl/5.10.1-fasrc05    #5.26.1-fasrc01

#module load gcc/7.1.0-fasrc01 mpich/3.2-fasrc03
#module load gcc/7.1.0-fasrc01 mvapich2/2.3b-fasrc02
#module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
module load gcc/7.1.0-fasrc01 netcdf/4.5.0-fasrc02
#module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02 netcdf-fortran/4.4.0-fasrc03
module load gcc/7.1.0-fasrc01 mvapich2/2.3b-fasrc02 netcdf-fortran/4.4.0-fasrc03

module load ncl_ncarg/6.4.0-fasrc01
module load matlab/R2018a-fasrc01
module load python/3.6.3-fasrc02
module load cdo/1.9.4-fasrc02
#module load gcc/7.1.0-fasrc01
module load ncview/2.1.7-fasrc02
#module load intel/17.0.4-fasrc01


module load git/2.17.0-fasrc01

# for GEOS-Chem
export FC=gfortran
export CC=gcc
export CXX=g++
export COMPILER=$FC

export NETCDF=/n/helmod/apps/centos7/Comp/gcc/7.1.0-fasrc01/netcdf/4.5.0-fasrc02
export GC_BIN=$NETCDF/bin
export GC_INCLUDE=$NETCDF/include
export GC_LIB=$NETCDF/lib

export NETCDF_FORTRAN=/n/helmod/apps/centos7/MPI/gcc/7.1.0-fasrc01/mvapich2/2.3b-fasrc02/netcdf-fortran/4.4.0-fasrc03
export GC_F_BIN=$NETCDF_FORTRAN_HOME/bin
export GC_F_INCLUDE=$NETCDF_FORTRAN_INCLUDE
export GC_F_LIB=$NETCDF_FORTRAN_LIB

# for WRF
export NETCDF_LIB=$NETCDF/lib:
export NETCDF_INC=$NETCDF/include:
export PATH=$NETCDF/bin:$PATH

export JASPERLIB=/n/helmod/apps/centos7/Comp/gcc/7.1.0-fasrc01/jasper/1.900.1-fasrc02/lib64:/n/helmod/apps/centos7/Core/libpng/1.6.25-fasrc01/lib:/n/helmod/apps/centos7/Core/zlib/1.2.8-fasrc09/lib
export JASPERINC=/n/helmod/apps/centos7/Comp/gcc/7.1.0-fasrc01/jasper/1.900.1-fasrc02/include:/n/helmod/apps/centos7/Core/libpng/1.6.25-fasrc01/include:/n/helmod/apps/centos7/Core/zlib/1.2.8-fasrc09/include

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export WRF_EM_CORE=1
export WRF_NMM_CORE=0
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
#export PATH=/n/helmod/apps/centos7/Comp/intel/17.0.4-fasrc01/netcdf/4.5.0-fasrc02/include


