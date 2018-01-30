export DEVROOT=`pwd`
export VERS=dbg
#export VERS=opt

export MESHSIM=/projects/tools/SimmetrixTest/12.0-171109dev
export PARALLEL=openmpi
export PATH=/usr/local/openmpi/1.10.6-gnu49-thread/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/openmpi/1.10.6-gnu49-thread/lib:$LD_LIBRARY_PATH


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
