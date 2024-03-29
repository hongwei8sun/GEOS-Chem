# .bashr

# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

ulimit -s unlimited
export OMP_STACKSIZE=500m
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

alias ls='ls --color=auto'
alias cp='cp -i'
alias rm='rm -i'
alias mv='mv -i'
alias python='python3.8'
alias sacct='sacct --format JobID,JobName,Partition,Elapsed,ReqMem,MaxRSS,AllocCPUs,TotalCPU,State'
alias gsu="git submodule update --init --recursive"

export LS_COLORS=$LS_COLORS:'rs=0:di=01;34:ln=01;36:mh=00:pi=47;33:so=01;35:do=01;35:bd=47;33;01:cd=47;33;01:or=47;31;01:su=37;41:sg=30;43:ca=30;41:tw=30;42:ow=34;42:st=37;44:ex=01;32:*.sh=01;92:*.m=01;32:*.tar=01;31:*.tgz=01;31:*.arj=01;31:*.taz=01;31:*.lzh=01;31:*.lzma=01;31:*.tlz=01;31:*.txz=01;31:*.zip=01;31:*.z=01;31:*.Z=01;31:*.dz=01;31:*.gz=01;31:*.lz=01;31:*.xz=01;31:*.bz2=01;31:*.bz=01;31:*.tbz=01;31:*.tbz2=01;31:*.tz=01;31:*.deb=01;31:*.rpm=01;31:*.jar=01;31:*.rar=01;31:*.ace=01;31:*.zoo=01;31:*.cpio=01;31:*.7z=01;31:*.rz=01;31:*.jpg=01;35:*.jpeg=01;35:*.gif=01;95:*.bmp=01;35:*.pbm=01;35:*.pgm=01;35:*.ppm=01;35:*.tga=01;35:*.xbm=01;35:*.xpm=01;35:*.tif=01;35:*.tiff=01;35:*.png=01;35:*.svg=01;35:*.svgz=01;35:*.mng=01;35:*.pcx=01;35:*.mov=01;35:*.mpg=01;35:*.mpeg=01;35:*.m2v=01;35:*.mkv=01;35:*.ogm=01;35:*.mp4=01;35:*.m4v=01;35:*.mp4v=01;35:*.vob=01;35:*.qt=01;35:*.nuv=01;35:*.wmv=01;35:*.asf=01;35:*.rm=01;35:*.rmvb=01;35:*.flc=01;35:*.avi=01;35:*.fli=01;35:*.flv=01;35:*.gl=01;35:*.dl=01;35:*.xcf=01;35:*.xwd=01;35:*.yuv=01;35:*.cgm=01;35:*.emf=01;35:*.axv=01;35:*.anx=01;35:*.ogv=01;35:*.ogx=01;35:*.aac=00;36:*.au=00;36:*.flac=00;36:*.mid=00;36:*.midi=00;36:*.mka=00;36:*.mp3=00;36:*.mpc=00;36:*.ogg=00;36:*.ra=00;36:*.wav=00;36:*.axa=00;36:*.oga=00;36:*.spx=00;36:*.xspf=00;36:*.nc=01;33:*.nc4=01;33:*.dat=01;33:*.txt=01;33:*.csv=01;33:'

module purge

# For python
module load Anaconda3/2019.10
source activate SUNenv
export PATH="$PATH:/n/home12/hongwei/.conda/envs/SUNenv/bin/python"

# For GEOS-Chem
module load intel/18.0.5-fasrc01 
module load openmpi/4.0.1-fasrc01 
module load netcdf/4.7.1-fasrc01
module load netcdf-fortran/4.5.2-fasrc01
module load cmake/3.17.3-fasrc01
module load gdb/8.1-fasrc01
module load emacs/26.1-fasrc01
module load git/2.17.0-fasrc01

export GC_BIN=$NETCDF_HOME/bin
export GC_INCLUDE=$NETCDF_HOME/include
export GC_LIB=$NETCDF_HOME/lib64

export GC_F_BIN=$NETCDF_FORTRAN_HOME/bin
export GC_F_INCLUDE=$NETCDF_FORTRAN_HOME/include
export GC_F_LIB=$NETCDF_FORTRAN_HOME/lib

export F90=$FC
export F77=$FC


module load ncl_ncarg/6.4.0-fasrc01
#module load cdo
