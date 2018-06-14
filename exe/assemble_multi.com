#!/bin/bash

mkdir -p me_files
filename='./data/blocks_multi.dat'
numfiles=1
numread=0
while read dir
do
    if [ ! -z "$dir" ]; then
	numread=$(($numread+1))
	if [ $numread == 1 ]; then
	    if [[ ! $dir == "update_en" ]]; then
		if [[ ! $dir == "update_en_var" ]]; then
		    if [[ ! $dir == "update_gr_en" ]]; then
			echo $dir > head.me
		    fi
		fi
	    fi
	    if [[ $dir == "update_en" ]]; then
		numread=0
		numfileread=$(($numfiles-1))
		filename=./me_files/SP_en_$numfileread.me
		while read en
		do		    
		    echo "well energy is "$en
		    en_well=$en
		done <  $filename
		filename='./data/blocks_multi.dat'
		echo "well energy is "$en_well
		newen=$(bc <<EOF 
                ($en_well - $en_ref) * 627.5
EOF
)
		echo "new energy is "$newen
		sed -ie "/ZeroEnergy/c\ZeroEnergy[kcal/mol] $newen" ./me_files/SP_$numfileread.me
	    fi

	    if [[ $dir == "update_gr_en" ]]; then
		numread=0
		numfileread=$(($numfiles-1))
		filename=./me_files/SP_en_$numfileread.me
		while read en
		do		    
		    echo "bimol energy is "$en
		    en_well=$en
		done <  $filename
		filename='./data/blocks_multi.dat'
		echo "bimol energy is "$en_well
		newen=$(bc <<EOF 
                ($en_well - $en_ref) * 627.5
EOF
)
		echo "new energy is "$newen
		sed -ie "/GroundEnergy/c\GroundEnergy[kcal/mol] $newen" ./me_files/SP_$numfileread.me
	    fi

	    if [[ $dir == "update_en_var" ]]; then
		numread=0
		numfileread=$(($numfiles-1))
		filename=./me_files/SP_en_$numfileread.me
		while read en
		do		    
		    en_var_ref=$en
		done <  $filename
		en_var_ref_kcal=$(bc <<EOF 
		    ($en_var_ref - $en_ref) * 627.5
EOF
)
		echo "variational reference energy is "$en_var_ref_kcal

		egrep Zero ./me_files/SP_$numfileread.me > sub_var.dat
		filename=./sub_var.dat
		numprog=0
		while read enline
		do		    
		    numprog=$(($numprog+1))
		    echo $enline > temp.log
		    sed -ie "s/ZeroEnergy\[kcal\/mol]/ /g" temp.log
		    filename=./temp.log
		    while read ensub
		    do		    
			entosub=$ensub
		    done <  $filename
		    filename=./sub_var.dat
		    if [ $numprog == 1 ]; then
			en_0_var=$entosub
		    fi
		    en_var_diff=$(bc <<EOF 
		    $en_var_ref_kcal+$entosub-$en_0_var
EOF
)
		    echo "en to sub is $entosub" 
		    sed -ie "s/$entosub/$en_var_diff/g" ./me_files/SP_$numfileread.me
		done <  $filename

		filename='./data/blocks_multi.dat'
	    fi
	fi
	if [ $numread == 2 ]; then
	    echo "numfiles is $numfiles"
	    echo "taking data from directory $dir"
	    cp -f $dir ./temp.me
	    sed -i '1d' temp.me
	    cat head.me temp.me > ./me_files/SP_$numfiles.me
	fi
	if [ $numread == 3 ]; then
	    cp -f $dir ./me_files/SP_en_$numfiles.me
	    echo "numfiles is $numfiles"
	    echo "taking data from directory $dir"
#	    echo "numread is $numread"
	    if [ $numfiles == 1 ]; then
		filename=./me_files/SP_en_$numfiles.me
		while read en
		do		    
		    echo "reference energy is "$en
		    en_ref=$en
		done <  $filename
		filename='./data/blocks_multi.dat'
		echo "reference energy is "$en_ref
	    fi
	    numfiles=$(($numfiles+1))
	    numread=0
	fi
    fi
done <  $filename

numfiles=$(($numfiles-1))

cat data/me_head.dat > me_ktp.inp

for npoints in $(seq "$numfiles" )
do
    cat me_ktp.inp  ./me_files/SP_$npoints.me > me_ktp_temp.me
    cp -f me_ktp_temp.me me_ktp.inp
done
echo  " End" >> me_ktp.inp

exit 0




