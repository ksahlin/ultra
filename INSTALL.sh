#!/bin/bash
set -e
if [ $# -ne 1 ]; then
        echo "Usage: `basename $0` <installation_path> "
        exit -1
fi
path=$1
mkdir -p $path

# create temporary installation folder
mkdir -p "temp_install_ultra"
cd "temp_install_ultra"

echo
echo "I'm in temporary folder" $PWD
echo 

# Create and activate a new environment called ultra
echo "SETTING UP CONDA ENVIRONMENT"
echo
conda create --yes -n ultra python=3.8 pip 
conda activate ultra

# Install uLTRA
echo
echo "INSTALLING ULTRA"
echo
pip install ultra-bioinformatics


# Install MEM finder slaMEM 
echo
echo "INSTALLING SLAMEM"
echo
installed_slamem=false
if ! command -v slaMEM &> /dev/null
then
    echo "slaMEM not found in path. Installing"
    git clone https://github.com/ksahlin/slaMEM.git
    cd slaMEM
    make 
    mv slaMEM $path
    echo
    echo "I have put slaMEM in:" $path " please make sure that this folder is in your path, or move slaMEM to your path"
    echo
    cd ..
    rm -rf slaMEM
    installed_slamem=true
fi

# Install aligner minimap2
echo
echo "INSTALLING MINIMAP2"
echo

installed_mm2=false
if ! command -v minimap2 &> /dev/null
then
    echo "minimap2 not found in path. Installing"
    git clone https://github.com/lh3/minimap2
    cd minimap2 && make
    mv minimap2 $path

    echo
    echo "I have put minimap2 in:" $path " please make sure that this folder is in your path, or move minimap2 to your path"
    echo
    cd ..
    rm -rf minimap2
    installed_mm2=true
fi

# Install Mummer

echo
echo "INSTALLING MUMMER"
echo

conda install --yes -c bioconda mummer

cd ..
rm -rf "temp_install_ultra"


echo
echo "FINISHED INSTALLATION"
echo
if [ "$installed_mm2" = true ] ; then
    echo "I have put minimap2 in:" $path " please make sure that this folder is in your path, or move minimap2 to your path"
fi

if [ "$installed_slamem" = true ] ; then
    echo "I have put slaMEM in:" $path " please make sure that this folder is in your path, or move slaMEM to your path"
fi

echo "Please activate the environment as 'conda activate ultra' before running uLTRA."
echo