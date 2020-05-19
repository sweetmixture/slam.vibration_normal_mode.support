#!/bin/bash

mkdir tmp
mkdir fdm
mkdir hess

mv geo.txt.*_"$1"* ./tmp
mv atom_* ./fdm
mv *hessian.dat ./hess

mkdir final

mv tmp  final
mv fdm  final
mv hess final

tar -czf final.tar.gz final
rm -rf final
