#! /bin/bash

########################################################################################
#       Parameter block begins                                                         #
########################################################################################

lambda=8.0        #list of `t Hooft dimensionless area
nconf=2005        #max no. of configurations in each  MPI stream to be generated. Set this to high value.
ntherm=5         #number of initial config (=nauto*nntherm trajectories) to skip for fermion measurments
nconfnotune=0

fm=1.0            #Wilon mass

iboolwrite=0      #iboolwrite=1 will write configs sequentially to GAUGE/, iboolwrite=0 will instead rewrite gaugeold files

hstart=1          #hot=1 or cold=0 start
nauto=5           #no. of traj between each written config ~ autocorrelation time


nmd=25            #initial no. of molecular dynamics steps
dt=0.04           #initial MD step size

nproc=1          #no. of MPI processes
########################################################################################
#       Parameter block end                                                            #
########################################################################################
                   

#the array a below has to be changed to (000 001  .... (value of nproc padded with appropritate number of 0s in front))

a=(000)

#array of ell's to be run in diff cores
glist=(${lambda})

mkdir -p GAUGE
mkdir -p EIG
mkdir -p MES

for iproc in `seq 0 $(($nproc-1))`
do
  if !(test -e inputfile${a[$iproc]})
  then
    echo 0 > inputfile${a[$iproc]}
#  else
#    mv gaugeold* GAUGE
  fi
done

for iproc in `seq 0 $(($nproc-1))`
do

istart=$(<inputfile${a[$iproc]})

echo ${glist[$iproc]} ${fm} 1.0 > inp${a[$iproc]}
echo 3500  1.E-4 3500 1.E-6 >>inp${a[$iproc]}

  if test -e gaugeold${a[$iproc]}
  then
    if test -e data${a[$iproc]}
    then                                 # input and data files exist: continuation run
      echo -1  -1 ${iboolwrite} >> inp${a[$iproc]}
    else                                 # input file exists but no data file: continuation run
      echo -1   1 ${iboolwrite} >> inp${a[$iproc]}
    fi
  else
    if [ "${istart}" -eq "0" ]                                   # input file does not exist: new run
    then
       echo ${hstart}  1 ${iboolwrite} >> inp${a[$iproc]}
    else
       if test -e data${a[$iproc]}
       then                                 # input and data files exist: continuation run
          echo -2  -1 ${iboolwrite} >> inp${a[$iproc]}
       else                                 # input file exists but no data file: continuation run
          echo -2   1 ${iboolwrite} >> inp${a[$iproc]}
       fi
    fi
  fi

echo ${nmd} ${nconf} ${nauto} 1 1 >> inp${a[$iproc]}
echo $(($iproc+1)) >> inp${a[$iproc]}
echo ${ntherm} >> inp${a[$iproc]}
echo ${nconfnotune} >> inp${a[$iproc]}
echo 0 >> inp${a[$iproc]}


echo " " >> outrun${a[$iproc]}
done

(../hmc)

exit
