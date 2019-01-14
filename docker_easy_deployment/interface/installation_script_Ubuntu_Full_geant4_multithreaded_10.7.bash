#!/bin/bash
set -e

#################
mkdir -p geant4 # directory were everything is built and installed
cd geant4
#############

########################## VARIABLES

##############  PROGRAMS' VERSIONS AND URLs : MAY CHANGE IN THE FUTURE
g4_version=10.7.p03
_g4_version=10.07.p03
folder_g4_version=Geant4-10.7.3
g4_url=("http://cern.ch/geant4-data/releases/geant4.${_g4_version}.tar.gz")

cmake_download_url=https://github.com/Kitware/CMake/releases/download/v3.14.3/cmake-3.14.3-Linux-x86_64.tar.gz

xerces_w_ver=xerces-c-3.2.2
xerces_arc=${xerces_w_ver}.tar.gz
xerces_url=("http://archive.apache.org/dist/xerces/c/3/sources/$xerces_arc")

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
base_dir=${PWD}
CMake_path=${base_dir}/cmake/bin/cmake

#
current_dir=$PWD

# Parameters
core_nb=`grep -c ^processor /proc/cpuinfo`



# Geant4
src_dir=$base_dir/source_geant4.${_g4_version}/
build_dir=$base_dir/geant4_build_${_g4_version}/
install_dir=$base_dir/geant4_install_${_g4_version}/
geant4_lib_dir=${install_dir}/lib/${folder_g4_version}/

# XERCES-C

xercesc_build_dir=($base_dir/build_xercesc_g4_${_g4_version}/)
xercesc_install_dir=($base_dir/install_xercesc_g4_${_g4_version}/)
xercesc_inc_dir=(${xercesc_install_dir}/include)
xercesc_lib_dir=(${xercesc_install_dir}/lib64/libxerces-c-3.2.so)

########## Creating folders

mkdir -p ${build_dir} # -p will create only if it does not exist yet
mkdir -p ${src_dir}
mkdir -p ${install_dir}

mkdir -p $xercesc_build_dir
mkdir -p $xercesc_install_dir

############# CHECK IF OS IS UBUNTU
echo "checking if OS is Ubuntu..."
# Checking if OS is Ubuntu
if [ -f /etc/os-release ]; then
    . /etc/os-release
    OS=$NAME
    VER=$VERSION_ID
fi

if [ ! "$OS" = "Ubuntu" ]; then
    echo "Error: OS is not Ubuntu. Script works only for Ubuntu. Aborting."
    exit 1
else
    echo "... OS is Ubuntu"
fi
#############

#########################################################################
############# CHECK IF DEPENDENCIES ARE SATISFIED, OTHERWISE INSTALL THEM

ubuntu_dependences_list=( "build-essential"
    "qtcreator"
    "cmake-qt-gui"
    "gcc"
    "g++"
    "gfortran"
    "zlib1g-dev"
    "libxerces-c-dev"
    "libx11-dev"
    "libexpat1-dev"
    "libxmu-dev"
    "libmotif-dev"
    "libboost-filesystem-dev"
    "libeigen3-dev"
    "libuuid1"
    "uuid-dev"
    "uuid-runtime"
)

entered_one_time=true

run_install()
{
    echo "Some missing dependencies were detected."
    ## Prompt the user
    if [ entered_one_time=true ]; then
        entered_one_time=false
        read -p "Do you have (root) sudo access ? [Y/n]. It is required to install missing dependencies: " answer
        ## Set the default value if no answer was given
        answer=${answer:N}
        if [[ $answer =~ [Nn] ]]; then
            echo "root access is required to install missing dependencies. Aborting."
            exit 1
        fi
    fi
    ## Prompt the user
    sudo add-apt-repository ppa:rock-core/qt4
    read -p "Do you want to install missing dependencies? [Y/n]: " answer
    ## Set the default value if no answer was given
    answer=${answer:Y}
    ## If the answer matches y or Y, install
    if [[ $answer =~ [Yy] ]]; then
        sudo apt-get install ${ubuntu_dependences_list[@]}
    else
        echo "Missing dependencies are required for proper compilation and installation. Aborting."
        exit 0
    fi
}


echo "checking dependencies..."

dpkg -s "${ubuntu_dependences_list[@]}" >/dev/null 2>&1 || run_install

echo "... dependencies are satisfied."

#########################################################################


#### XERCES-C (to be able to use GDML files)

## download xerces-c (for GDML)

wget $xerces_url
tar zxf $base_dir/$xerces_arc
rm -rf $xerces_arc

xerces_src=$base_dir/$xerces_w_ver

## compile and install xerces-c

cd $xercesc_build_dir

echo "build of xerces-c: Attempt to execute CMake..."

rm -rf CMakeCache.txt

$CMake_path \
-DCMAKE_INSTALL_PREFIX=${xercesc_install_dir} \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_INSTALL_LIBDIR=lib64 \
$xerces_src
echo "... done"

echo "Attempt to compile and install xerces-c"

G4VERBOSE=1 make -j${core_nb}
make install

cd $base_dir
echo "... done"

#### GEANT4

## download Geant4

rm -rf ${src_dir}
wget $g4_url
tar zxf geant4.${_g4_version}.tar.gz
mv geant4.${_g4_version} ${src_dir}
rm -rf geant4.${_g4_version}.tar.gz

## compile and install Geant4

cd ${build_dir}
rm -rf CMakeCache.txt

echo "build_geant4: Attempt to execute CMake"

$CMake_path \
-DCMAKE_PREFIX_PATH=${xercesc_install_dir} \
-DCMAKE_INSTALL_PREFIX=${install_dir} \
-DCMAKE_BUILD_TYPE=Release \
-DGEANT4_BUILD_MULTITHREADED=ON \
-DGEANT4_INSTALL_DATA=ON \
-DGEANT4_USE_GDML=ON \
-DGEANT4_USE_G3TOG4=ON \
-DGEANT4_USE_QT=ON \
-DGEANT4_FORCE_QT4=OFF \
-DGEANT4_USE_XM=ON \
-DGEANT4_USE_OPENGL_X11=ON \
-DGEANT4_USE_INVENTOR=OFF \
-DGEANT4_USE_RAYTRACER_X11=ON \
-DGEANT4_USE_SYSTEM_CLHEP=OFF \
-DGEANT4_USE_SYSTEM_EXPAT=OFF \
-DGEANT4_USE_SYSTEM_ZLIB=OFF \
-DCMAKE_INSTALL_LIBDIR=lib \
-DXERCESC_INCLUDE_DIR=${xercesc_inc_dir} \
-DXERCESC_LIBRARY=${xercesc_lib_dir} \
../source_geant4.${_g4_version}/

echo "... Done"

echo "Attempt to compile and install Geant4"

G4VERBOSE=1 make -j${core_nb}

make install

cd $base_dir
echo "... Done"

#########################################################################
#########################################################################
#### set environement variables into '~/.bashrc'

echo "Attempt to setup up environement variables..."

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

# clean environement that was previously set by this script
first_line=`grep -n "## --> Added by Geant4 installation script" ~/.bashrc | awk -F  ":" '{print $1}'`
echo $first_line
last_line=`grep -n "## <-- Added by Geant4 installation script" ~/.bashrc | awk -F  ":" '{print $1}'`
echo $last_line

re='^[0-9]+$'
if [[ $first_line =~ $re ]] ; then # if $first_line is a number (i.e. it was found)
    if [[ $last_line =~ $re ]] ; then # if $last_line is a number (i.e. it was found)
        sed -i.bak "${first_line},${last_line}d" ~/.bashrc # delete text in .bashrc from first-line to last-line
    fi
fi

#
echo "## --> Added by Geant4 installation script" >> ~/.bashrc

set_environement() {
    
    cd ${base_dir}
    
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

# Geant4 + data
set_environement "source ${install_dir}/bin/geant4.sh"

# xerces-c
set_environement "export C_INCLUDE_PATH=\$C_INCLUDE_PATH:${xercesc_install_dir}/include/"
set_environement "export CPLUS_INCLUDE_PATH=\$CPLUS_INCLUDE_PATH:${xercesc_install_dir}/include/"
set_environement "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:${xercesc_install_dir}/lib64/"
set_environement "export LIBRARY_PATH=\$LIBRARY_PATH:${xercesc_install_dir}/lib64/"
set_environement "export PATH=\$PATH:${xercesc_install_dir}/include/"

echo " " >> ~/.bashrc
echo "## <-- Added by Geant4 installation script" >> ~/.bashrc
echo "... Done"
echo -e "${RED}Please excecute command < ${GREEN}source ~/.bashrc${RED} > or re-open a terminal for the system to be able to find the databases and libraries.${NC}"


