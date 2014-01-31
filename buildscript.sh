#!/bin/bash -l

function usage {
    echo "bla bla bla"
    exit 1
}

BASEDIR=$(pwd)

TARGET=cpu
BUILDDIR="${BASEDIR}/build"
STELLADIR="stella"

# Getting machine
machine=$(hostname | sed -rn 's/^([a-z]+).*/\1/p')
if [[ $machine == "ppopcode" ]]
then
    machine=opcode
fi

# Parsing options
while getopts "t:b:s:" o
do
    case "${o}" in
        t) # Target
            if [[ ${OPTARG} == "gpu" ]]
            then
                TARGET=gpu
            elif [[ ${OPTARG} == "cpu" ]]
            then
                TARGET=cpu
            else
                usage
            fi
            ;;

        b) # Build directory
            BUILDDIR=$(readlink --canonicalize "${OPTARG}")
            ;;
        
        s) # Stencil library
            STELLADIR=$(readlink --canonicalize "${OPTARG}")
            ;;
    esac
done

CMAKEARGS=" \
    -DCMAKE_BUILD_TYPE=Release \
"

# Programming environment
loadedPrgEnv=$(module list 2>&1 | sed -rn 's|.*(PrgEnv-[a-z]+).*|\1|p')
module switch $loadedPrgEnv PrgEnv-gnu
export CC=cc
export CXX=CC
module load cmake

# Boost
if [[ $machine == "todi" ]]
then
    module load boost/1.51.0
elif [[ $machine == "daint" ]]
then
    export BOOST_ROOT=/apps/daint/boost/1.54.0/gnu_473
fi

# CUDA
if [[ $TARGET == "gpu" ]]
then
    CMAKEARGS+=" \
        -DCUDA_BACKEND=ON \
        -DCUDA_NVCC_FLAGS=-I${MPICH_DIR}/include \
    "
    module load cudatoolkit
else
    CMAKEARGS+="-DCUDA_BACKEND=OFF  "
fi

# STELLA
CMAKEARGS+="-DSTELLA_PATH=${STELLADIR}  "

# Enter directory
mkdir -p $BUILDDIR
cd $BUILDDIR

# Call cmake
if [[ ! -a CMakeCache.txt ]]
then
    echo "Calling cmake $BASEDIR $CMAKEARGS" 
    cmake "${BASEDIR}" $CMAKEARGS
fi

if [[ $? -gt 0 ]]
then
    echo "Error configuring. Aborting"
    exit 2
fi

# Build
make -j8
if [[ $? -gt 0 ]]
then
    echo "Error building. Aborting"
    exit 3
fi

