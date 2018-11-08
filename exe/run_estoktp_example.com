#!/bin/sh

SSH=/usr/bin/ssh

HOST=$1
#NODESLOW=$2
#NODESHIGH=$3

DIRECTORY='./data'

if [ ! -d "$DIRECTORY" ] 
then
    echo " "
    echo " it seems the subdirectory data is missing"
    echo " check and run again"
    exit 
fi


# if [ "$#" != "3" ] 
if [ "$#" != "1" ] 
then
#    echo " usage  ./run_estoktp.com  node #procLowLev #procHighLev"
    echo " usage  ./run_estoktp.com  node "
fi

#if [ "$#" != "3" ] 
#then
#    echo " input file in ./data/estoktp.dat"
#fi

#if [ "$#" != "3" ] 
if [ "$#" != "1" ] 
then
    exit 
fi

echo " "
#echo "nodes low level = $NODESLOW"
#echo "nodes high level = $NODESHIGH"
echo "host = $HOST"
echo "$HOST" > host.dat

#sed -ie 's/NODESLOW/'$NODESLOW'/g' data/estoktp.dat
#sed -ie 's/NODESHIGH/'$NODESHIGH'/g' data/estoktp.dat



# bash command
#  exec $SSH -n $HOST "cd `pwd`; export PATH=$PATH:~/bin; estoktp.x >& estoktp.log & "

# csh tcsh comman
  exec $SSH -n $HOST "cd `pwd`; setenv PATH $PATH:~/bin; estoktp.x >& estoktp.log & "




