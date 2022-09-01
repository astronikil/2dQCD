#! /bin/bash

# commit all files #

echo "Check if the entries in header file is correct:"

echo ">>>>>>>> BEGIN HEADER <<<<<<<<<<<<<<<"
cat ../headers/size.h
echo ">>>>>>>>   END HEADER <<<<<<<<<<<<<<<"


../hmc/comp.x

echo 'Compiled HMC'

../ritz_rajamani/comp.x
echo 'Compiled eigen'

../c2pt/comp.x
echo 'Compiled correl'



