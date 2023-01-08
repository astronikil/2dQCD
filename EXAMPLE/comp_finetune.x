#! /bin/bash

# commit all files #

rm -f hmc 

echo " "
echo " "
echo "Check if the entries in header file is correct:"
echo " "
echo " "

echo " "
echo " "
echo ">>>>>>>> BEGIN HEADER <<<<<<<<<<<<<<<"
cat ../headers/size.h
echo ">>>>>>>>   END HEADER <<<<<<<<<<<<<<<"
echo " "
echo " "

../hmc/comp_finetune.x


echo "Compile HMC done"


