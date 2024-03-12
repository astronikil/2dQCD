#! /bin/bash

rm -f eig.dat

./suneig
cat fort.55 >> eig.dat

