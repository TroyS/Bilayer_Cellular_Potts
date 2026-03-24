#!/bin/bash


fnew=$1
fin=$2
fpara=$3
fout=$4

g++ -o Bilaybisweep BilayerCPbisweep.cpp

cp Bilaybisweep $fnew
cp $fpara $fnew
cp cdat2vtkfull $fnew  
cd $fnew

./Bilaybisweep $fin $fpara $fout

