export DEVROOT=`pwd`
#export PARALLEL=mpich2
export PARALLEL=bgmpi
export ARCHOS=ibm-bg
export VERS=opt
#export MESHSIM=/usr/local/simmetrix/simmodsuite/9.0-130610
export MESHSIM=/home/chitak/meshSim/9.0-130517/
#export MESHSIM=/usr/local/simmetrix/simmodsuite/7.2-120626
#export MESHSIM=/users/chitak/meshSim/main-130220
##source env.sh
#export MPIHOME=/soft/compilers/wrappers/xl
export CC=mpixlc
export CXX=mpixlcxx
export FC=mpixlf90
#export CC=mpicc
#export CXX=mpicxx
#export FC=gfortran

#isclean='distclean'
compile="make SIM=1 MODELER=geomsim NODEP=1 NOSHARED=1 $isclean"
setup="make SIM=1 MODELER=geomsim NODEP=1 NOSHARED=1 setup"

dest_path=$DEVROOT/phasta/phastaIO/phastaIO
cd $dest_path
$setup
$compile

dest_path=$DEVROOT/phasta/phShape/phShape
cd $dest_path
$setup
$compile

dest_path=$DEVROOT/phasta/phUtil/LU/LU
cd $dest_path
$compile

dest_path=$DEVROOT/phParAdaptRed/phParAdapt
cd $dest_path
$compile
