#!/bin/bash

########
# Utility script for the purpose of applying possible displacements of atoms within imaginary vibration frequencies
#
# Simply ...
# $ bash apply_displs.sh [slam geometry file (geo.txt)]  [mode number]
#
# Make sure both
# 'vib_displ.out' : displacement log - standard output from SLAM normal model calculation
# 'geo.txt'       : standard SLAM geometry input file (target to apply displacements)
# files are in the same directory
########

TARGET_GEO="geo.txt"
DISPLS_INFO="vib_displ.out"
TARGET_MODE=$1

DEFAULT_W="1.0"

rm -rf geo.txt.displ
touch geo.txt.displ
########
# CHECK IF FILE EXISTS !
########

if [ -z "$TARGET_GEO" ]; then
    echo "Error, geo.txt does not exists, please do double check!"
    exit 1
fi
if [ -z "$DISPLS_INFO" ]; then
    echo "Error, vib_displ.out does not exists, please do double check!"
    exit 1
fi
if [ -z "$TARGET_MODE" ]; then
    echo "Error, which mode you want to apply? second argument is not specified!"
    exit 1
fi

# TRY catch corresponding mode
tar="normal_mode_# ""$TARGET_MODE"
IF_TC_SUCCESS=$( grep "$tar" $DISPLS_INFO )
if [ -z "IF_TC_SUCCESS" ]; then
    echo "Error, specified mode" "$TARGET_MODE" "does not exist!"
    exit 1
else
    NOA=$( grep -B 1 "$tar" $DISPLS_INFO | head -1 )
    grep -A "$NOA" "$tar" $DISPLS_INFO | tail -n $NOA > target_displs_tmp       # tempolar target displs 
fi  

declare -a displs

for(( i=0; i<"$NOA"; i++ )); do
    ln=$(echo "$i+1" | bc)
    line=$(sed -n "$ln"p "target_displs_tmp")
    spl=( $line )
    for (( j=0; j<4; j++ )); do
        displs[$i*4+$j]=${spl[$j]}
    done
done


#for(( i=0; i<"$NOA"; i++ )); do
#    echo ${displs[$i*4+0]} ${displs[$i*4+1]} ${displs[$i*4+2]} ${displs[$i*4+3]}
#done

# SAVE DISPLS IN ARR 'displs'
head -2 geo.txt > geo.txt.displ

i=0     # WORKING INDICES
j=0
while IFS= read -r line
do
    i=$((i+1))
    if [ $i -gt 2 ]; then
        spl=( $line )

        re='^[+-]?[0-9]+([.][0-9]+)?$'  # check if a variable is a number

        if [ ${spl[1]} == 'c' ]; then

            x_org=${spl[2]}
            x_del=${displs[$j*4+1]}
            y_org=${spl[3]}
            y_del=${displs[$j*4+2]}
            z_org=${spl[4]}
            z_del=${displs[$j*4+3]}     # if target is classical core part + add displacement
        
            x=$(echo "$x_org + $x_del" | bc -l )
            y=$(echo "$y_org + $y_del" | bc -l )
            z=$(echo "$z_org + $z_del" | bc -l )

            j=$((j+1))
            printf "%.3s%2.1s%9.6f%12.6f%12.6f\n" ${spl[0]} ${spl[1]} $x $y $z >> geo.txt.displ # SAVE CLASSICAL ION

        elif [ ${spl[1]} == 's' ]; then
            echo $line >> geo.txt.displ # if target is 'shell' + just add

        elif [[ ${spl[1]} =~ $re ]]; then       # if target LP ion ... 
            
            x_org=${spl[1]}
            x_del=${displs[$j*4+1]}
            y_org=${spl[2]}
            y_del=${displs[$j*4+2]}
            z_org=${spl[3]}
            z_del=${displs[$j*4+3]}
        
            x=$(echo "$x_org + $x_del" | bc -l )
            y=$(echo "$y_org + $y_del" | bc -l )
            z=$(echo "$z_org + $z_del" | bc -l )

            j=$((j+1))
            printf "%.3s%9.6f%12.6f%12.6f\n" ${spl[0]} $x $y $z >> geo.txt.displ        # SAVE LP ION
        fi
    fi
done < "geo.txt"

# clean temporals
rm -rf target_displs_tmp
