#! /bin/bash

rm -r fort.* *000 GAUGE


fm=1.0            #Wilon mass

iboolwrite=0   #iboolwrite=1 will write configs sequentially to GAUGE/, iboolwrite=0 will instead rewrite gaugeold files

stringlen=0

hstart=0          #hot=1 or cold=0 start
nauto=5           #no. of traj between each written config

g1=4.0           #list of 4 different ell's with 32 streams each
g2=4.0
g3=4.0
g4=4.0

itune=0           #to get a detailed report on time taken in trajectories
nmd=39             #initial no. of molecular dynamics steps
ntherm=0          #put 0; doesn't matter
dt=0.025            #initial MD step size

nproc=1        #no. of MPI processes
nconf=200        #max no. of configurations in each  MPI stream to be generated. Set this to high value.
                   

#the array a below has to be changed to (000 001  .... (value of nproc padded with appropritate number of 0s in front))

a=(000)



#array of ell's to be run in diff cores
glist=(${g1})

mkdir GAUGE

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

echo ${glist[$iproc]} ${fm} > inp${a[$iproc]}
echo 3500  1.E-5 3500 1.E-7 >>inp${a[$iproc]}

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

echo " " >> outrun${a[$iproc]}
done

(../su2)

exit

