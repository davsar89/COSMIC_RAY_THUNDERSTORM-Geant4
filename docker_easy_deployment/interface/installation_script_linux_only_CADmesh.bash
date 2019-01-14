#!/bin/bash
set -e

#################
mkdir -p geant4 # directory were everything is built and installed
cd geant4
base_dir=$PWD
#############
source /geant4/geant4_install_10.07.p03/bin/geant4.sh || true
git config --global http.sslverify "false" || true
########################## VARIABLES

##############  PROGRAMS' VERSIONS AND URLs : MAY CHANGE IN THE FUTURE

cmake_download_url=https://github.com/Kitware/CMake/releases/download/v3.14.3/cmake-3.14.3-Linux-x86_64.tar.gz

casmesh_w_ver=1.1
casmesh_arc=v${casmesh_w_ver}.tar.gz
casmesh_url=https://github.com/davsar89/CADMesh.git

####################################################

# getting CMake
rm -rf cmake
rm -rf cmake-3.14.3-Linux-x86_64
rm -rf cmake-3.14.3-Linux-x86_64.tar.gz
wget ${cmake_download_url}
tar zxf cmake-3.14.3-Linux-x86_64.tar.gz
mv cmake-3.14.3-Linux-x86_64 cmake
rm -rf cmake-3.14.3-Linux-x86_64.tar.gz

# CMake command
CMake_path=${base_dir}/cmake/bin/cmake

#
current_dir=$PWD

# Parameters
core_nb=`grep -c ^processor /proc/cpuinfo`

# CADMESH

casmesh_build_dir=($base_dir/build_cadmesh_g4_${_g4_version}/)
casmesh_install_dir=($base_dir/install_cadmesh_g4_${_g4_version}/)

rm -rf $casmesh_build_dir
rm -rf $casmesh_install_dir

########## Creating folders

mkdir -p $casmesh_build_dir
mkdir -p $casmesh_install_dir

#### CADMESH
# CADMESH is a CAD file interface for GEANT4, made by Poole, C. M. et al.
# See https://github.com/christopherpoole/CADMesh

## download CADMESH

rm -rf CADMesh

git clone $casmesh_url
casmesh_src=$base_dir/CADMesh

## compile and install CADMESH

cd $casmesh_build_dir

echo "build of CADMESH: Attempt to execute CMake..."

rm -rf CMakeCache.txt

$CMake_path \
-DCMAKE_INSTALL_PREFIX=${casmesh_install_dir} \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_INSTALL_LIBDIR=lib \
-DGeant4_DIR=$geant4_lib_dir \
$casmesh_src

echo "... done"

echo "Attempt to compile and install CADMESH"

G4VERBOSE=1 make -j${core_nb}
make install

cd $base_dir

echo "... done"

#########################################################################
#########################################################################
#### set environement variables into '~/.bashrc'

echo "Attempt to setup up environement variables..."

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

# clean environement that was previously set by this script
first_line=`grep -n "## --> Added by CADmesh installation script" ~/.bashrc | awk -F  ":" '{print $1}'`
echo $first_line
last_line=`grep -n "## <-- Added by CADmesh installation script" ~/.bashrc | awk -F  ":" '{print $1}'`
echo $last_line

re='^[0-9]+$'
if [[ $first_line =~ $re ]] ; then # if $first_line is a number (i.e. it was found)
    if [[ $last_line =~ $re ]] ; then # if $last_line is a number (i.e. it was found)
        sed -i.bak "${first_line},${last_line}d" ~/.bashrc # delete text in .bashrc from first-line to last-line
    fi
fi

#

echo "## --> Added by hdf5 zlib matio installation script" >> ~/.bashrc

set_environement() {
    
    cd $base_dir
    
    if grep -Fxq "$1" ~/.bashrc
    then
        echo -e "${GREEN}< source $1 > already set up in ~/.bashrc.${NC}"
    else
        echo "    " >> ~/.bashrc
        echo $1 >> ~/.bashrc
        echo "______"
        echo -e "${GREEN}added ${RED}$1${GREEN} to ${RED}~/.bashrc${GREEN} file.${NC}"
    fi
}

# CADMesh
set_environement "export cadmesh_DIR=${casmesh_install_dir}/lib/cmake/cadmesh-1.1.0/"
set_environement "export C_INCLUDE_PATH=\$C_INCLUDE_PATH:${casmesh_install_dir}/include/"
set_environement "export CPLUS_INCLUDE_PATH=\$CPLUS_INCLUDE_PATH:${casmesh_install_dir}/include/"
set_environement "export PATH=\$PATH:${casmesh_install_dir}/include/"
set_environement "export LIBRARY_PATH=\$LIBRARY_PATH:${casmesh_install_dir}/lib/"
set_environement "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:${casmesh_install_dir}/lib/"

echo "## <-- Added by CADmesh installation script" >> ~/.bashrc

echo "... Done"
echo -e "${RED}Please excecute command < ${GREEN}source ~/.bashrc${RED} > or re-open a terminal for the system to be able to find the databases and libraries.${NC}"



