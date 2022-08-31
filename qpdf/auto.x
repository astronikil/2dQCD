#! /bin/bash 

awk '{if($4=="S0V3S0"){print}}' fort.37  > foo.dat
awk '{print $1, $2}' foo.dat > foo3
awk '{if($4=="x23=" && $6=="S0"){print}}' fort.36  > g1
awk '{if($4=="x31=" && $6=="S0"){print}}' fort.36  > g2
awk '{if($4=="x12=" && $6=="S0"){print}}' fort.36  > g3
awk '{print $2}' g1 > foo1
awk '{print $2}' g2 > foo2
paste foo3 foo1 foo2  > out

