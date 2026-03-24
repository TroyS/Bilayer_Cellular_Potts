#!/bin/bash

fin=$1

fout=$2

i1=$3
i2=$4
i3=$5

g++ -o startBilay BilayerCPinit.cpp

./startBilay $fin $fout

./startBilay $fout $i1 $i2 $i3

mkdir $fout 

mv $fout* $fout
