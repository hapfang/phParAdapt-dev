export DEVROOT=`pwd`
#export VERS=dbg
export VERS=opt

#BGQ
#export PARALLEL=bgmpi
#export ARCHOS=ibm-bg
#export MESHSIM=/home/chitak/meshSim/9.0-130517/

#export CC=mpixlc
#export CXX=mpixlcxx
#export FC=mpixlf90

# UCB
export PARALLEL=openmpi
export ARCHOS=x86_64_linux
export MESHSIM=/users/mrasquin/develop/Meshing/simmodsuite-9.0-140927/9.0-140927
+simmodsuite-9.0-150808
export CC=mpicc
export CXX=mpicxx
export FC=gfortran

#export CXXFLAGS='-fsanitize=address -fno-omit-frame-pointer -g'
#export CFLAGS='-fsanitize=address -fno-omit-frame-pointer -g'
#export LDFLAGS='-fsanitize=address -fno-omit-frame-pointer -g'

#isclean='distclean'
#compile="make SIM=1 MODELER=geomsim NODEP=1 NOSHARED=1 $isclean"
#setup="make SIM=1 MODELER=geomsim NODEP=1 NOSHARED=1 setup"

compile="make SIM=1 MODELER=parasolid NODEP=1 NOSHARED=1 $isclean"
setup="make SIM=1 MODELER=parasolid NODEP=1 NOSHARED=1 setup"

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

dest_path=$DEVROOT/phParAdapt-Sim/phParAdapt
cd $dest_path
$compile
