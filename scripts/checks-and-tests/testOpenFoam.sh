#! /bin/bash

cd ..
rm -rf Yade-OpenFOAM-coupling

if [ -z "$WM_PROJECT_VERSION" ]; then
    if [ -f "/root/OpenFOAM/OpenFOAM-v1906/etc/bashrc" ]; then
        bashrcPath=/root/OpenFOAM/OpenFOAM-v1906/etc/bashrc # compiled version (ubuntu18.04)
    elif [ -f "/usr/lib/openfoam/openfoam2312/etc/bashrc" ]; then
        bashrcPath=/usr/lib/openfoam/openfoam2312/etc/bashrc # precompiled package (ubuntu22.04)
    else #assume OFOAM6, use older coupling code
        bashrcPath=/root/OpenFOAM/OpenFOAM-6/etc/bashrc
    fi
    source $bashrcPath
fi



cp -rf trunk/pkg/openfoam/coupling Yade-OpenFOAM-coupling
cd Yade-OpenFOAM-coupling
python3 setup.py

cd ../trunk/examples/openfoam/example_icoFoamYade
echo `pwd`
blockMesh
cp -r 0_org 0
decomposePar
mkdir spheres
#Create a symbolic link to Yade
# ln -s ../../../install/bin/yade-ci yadeimport.py


if which icoFoamYade &> /dev/null; then
    echo 'Command exists.'
else
    echo 'Command does not exist.'
fi

mpirun --allow-run-as-root -n 4 ../../../install/bin/yade-ci scriptMPI.py

echo -e "******************************************\n*** icoFoamYade test finished ***\n******************************************\n"


cd ../example_pimpleFoamYade
echo `pwd`
blockMesh
cp -r 0_org 0
decomposePar
mkdir spheres
# #Create a symbolic link to Yade
# ln -s ../../../install/bin/yade-ci yadeimport.py

if which pimpleFoamYade&> /dev/null; then
    echo 'Command exists.'
else
    echo 'Command does not exist.'
fi

mpirun --allow-run-as-root -n 4 ../../../install/bin/yade-ci scriptMPI.py

echo -e "******************************************\n*** pimpleFoamYade test finished ***\n******************************************\n"
