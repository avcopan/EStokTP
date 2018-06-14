#!/bin/bash

filename='./data/rea_pro_1.dat'
numreac=0
while IFS=' ' read -r name  reac
do
    if [ ! -z "$name" ]; then
	numreac=$(($numreac+1))
	declare name$numreac=${name}
	declare reac$numreac=${reac}
#	name_1=${name}
#	reac_1=${reac}
    fi
done <  $filename

#exit
#name1='ch3oh'
#reac1='/home/carlo/qm/combustion/database/C1/ch3oh'

cp -f ${reac1}/data/tsr1.dat .
cp -f data/rea_pro_2.dat .
#exit
#egrep nsites tsr1.dat > tmp.dat
filename='./tsr1.dat'
numsite=0
while IFS=' ' read -r junk
do
    if [ ! -z "$junk" ]; then
	if [[ $junk == "site" ]]; then
	    numsite=$(($numsite+1))
	    echo "site "$numsite > tsr_$numsite.dat
	fi
	if [[ ! $junk == "site" ]]; then
	    echo $junk  > tmp1.dat
	    cat tsr_$numsite.dat tmp1.dat  > tmp2.dat
	    cp -f tmp2.dat tsr_$numsite.dat
	fi
    fi
done <  $filename

rm -f tmp1.dat
rm -f tmp2.dat

#echo $prod1 > file1.dat

#exit

#prod1='/home/carlo/qm/combustion/ch3cooh/ch3coo'
filename='./rea_pro_2.dat'
numfiles=1
numread=0

for npoint in $(seq "$numsite" )
do
    while IFS=' ' read -r name2 dir_r2 dir_p2
    do
	if [ ! -z "$name2" ]; then
	    numread=$(($numread+1))
#	    echo  $dir_r2 > file1.dat
#	    echo  $dir_p2 > file2.dat
	    namedir=${name1}_${name2}_site${npoint}
#	    echo $namedir > file3.dat
#	    exit
	    mkdir -p $namedir
	    cd $namedir
	    mkdir -p data
	    mkdir -p output
	    mkdir -p geoms
	    mkdir -p irc_files
	    mkdir -p me_files
	    mkdir -p me_blocks
	    mkdir -p hl_logs
	    mkdir -p hr_geoms
	    mkdir -p md_tunn
#	    which get_files_nc.com > file4.dat
	    prod1=$(egrep prod1 ../tsr_$npoint.dat)
	    prod1="${prod1//prod1=}" 
#	    prod1=egrep prod1 ../tsr_$npoint.dat
	    get_files_nc.com reac1 reac1 $reac1
	    get_files_nc.com reac2 reac1 $dir_r2
	    get_files_nc.com prod1 reac1 $prod1
	    get_files_nc.com prod2 reac1 $dir_p2

	    sed -ie 's/Well REACS/Bimolecular REACS/g' me_files/reac1_ge.me
	    sed -ie 's/Species/Fragment REACT1/g' me_files/reac1_ge.me
	    sed -ie 's/Species/Fragment REACT2/g' me_files/reac2_ge.me
	    sed -ie '/End/d' me_files/reac1_fr.me
            cp -f ./me_files/reac1_fr.me .
            echo ' End' > temp.log
            echo '!*********' >> temp.log
            cat reac1_fr.me temp.log > ./me_files/reac1_fr.me

	    sed -ie '/End/d' me_files/reac2_fr.me
#	    sed -ie 's/End/ /g' me_files/reac2_fr.me
	    sed -ie 's/GroundEnergy\[kcal\/mol\] 0.0/ /g' me_files/reac2_fr.me
            cp -f ./me_files/reac2_fr.me .
            echo ' End' > temp.log
            echo '!*********' >> temp.log
            echo 'GroundEnergy[kcal/mol] 0.0' >> temp.log
            echo 'End' >> temp.log
            cat reac2_fr.me temp.log > ./me_files/reac2_fr.me

	    sed -ie 's/Well REACS/Bimolecular PRODS/g' me_files/prod1_ge.me
	    sed -ie 's/Species/Fragment PROD1/g' me_files/prod1_ge.me
	    sed -ie 's/Species/Fragment PROD2/g' me_files/prod2_ge.me
	    sed -ie '/End/d' me_files/prod1_fr.me
            cp -f ./me_files/prod1_fr.me .
            echo ' End' > temp.log
            echo '!*********' >> temp.log
            cat prod1_fr.me temp.log > ./me_files/prod1_fr.me
	    sed -ie '/End/d' me_files/prod2_fr.me
	    sed -ie 's/GroundEnergy\[kcal\/mol\] $proden/ /g' me_files/prod2_fr.me
            cp -f ./me_files/prod2_fr.me .
            echo ' End' > temp.log
            echo 'GroundEnergy[kcal/mol] $proden' >> temp.log
            echo 'End' >> temp.log
            echo '!*********' >> temp.log
            cat prod2_fr.me temp.log > ./me_files/prod2_fr.me

	    cp -f ../tsr_$npoint.dat ./data/tsr1.dat
	    cp -f $dir_r2/data/tsr2.dat ./data
	    cp -f ../data/estoktp.dat ./data
	    cp -f ../data/theory.dat ./data
	    cp -f ../data/hl_molpro.dat ./data
	    cp -f ../data/hl_ts_molpro.dat ./data
	    cp -f ../data/me_head.dat ./data
	    cd data
	    tsbuild.exe
	    cd ../..
	    rm -f  tsr_$npoint.da

	fi
    done <  $filename

done

rm -f  rea_pro_2.dat
rm -f  tsr1.dat
