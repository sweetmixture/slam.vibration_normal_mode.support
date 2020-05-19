#!/bin/bash

# BEWARE ... BC CANNOT HANDLE 'e-10' LIKE NUMBERS

###############################################################################################
# Last time editied : 29 - 03 - 2020
# Written by        : Woongkyu Jee
# Affiliation       : University College London, Dept of Chemistry.
# Contact	    : woong.jee.16@ucl.ac.uk / wldndrb1@gmail.com
#
# Description       : The main purpose of this script is for generating input files 'geo.txt'
# for numerical calculations of Hessian matrix of SLAM optimised clusters (molecules). As the
# numerical treatment, finite differnce method has introduced, which is,
#
#				F( x + dx ) - F( x - dx )
#		H = 	----------------------------------------
#					  2 dx
#
# For the variation 'dx' for calculation, users must think about trade-off versus the accuracy
# of Hessian figures and the following numerical noises. For running this script, in the follo
# -ing part of the comments, the default environmetal variables are set. Users have to adapt
# those, paths or file names etc., for your working or computer enviroments.
# For running calculations, keep potential (sp_cluster_species.txt) and basis set directory
# in the same working directory with 'geo.txt'. 
###############################################################################################
printf "%s\n" "######################################################################################"
printf "%s\n" "--------------------------------------------------------------------------------------"
printf "%s\n" " *supportive script for calculating slam normal modes" 
printf "%s\n" " compatibilty : SLAM V 2.0 with RIM calculations"
printf "%s\n" " last edited  : 29 - 03 -2020"
printf "%s\n" "--------------------------------------------------------------------------------------"

# Default env - variables.
TAR="geo.txt"								# Takes standard slam input geometry file 'geo.txt'
DELTA=$1								# Takes "dx" value, if it is not specified. then 								

printf "\n"
if [ -z $DELTA ]; then							# falling back to the default value '0.0025'
	DELTA="0.0025"							# default delta value ... 
	printf "%2s %12.6f\n" " Default 'dx' value is used :" $DELTA	#
else									#
	printf "%2s %12.6f\n" " Requested 'dx' value is used :" $DELTA	#
fi									#

# SLAM bin path & excutable ... user customisation available
SLAM_BINDIR="/home/uccawkj/update_slam_backup/bin"			# SLAM BIN DIR
SLAM_BIN="slam_v_2.0_opti.mpi.x"					# SLAM BIN PATH
PREFIX="mpirun -np"							# Parallel process arguments
NCORE="4"								# Number of cores requested for the calculation (single point)
EXE="$PREFIX $NCORE $SLAM_BINDIR/$SLAM_BIN"
DIGPY="/home/uccawkj/update_slam_backup/bin/slam_vib_package_v_2.0/slam_vib.diagonaliser.py"
CLEAN="/home/uccawkj/update_slam_backup/bin/slam_vib_package_v_2.0/clean.sh"
# ENV Variables check
if [ -z $SLAM_BINDIR ]; then									# Check SLAM_BINDIR correct
	echo "SLAM_BINDIR is not set ... check the coressponding variable !"
	exit 1
fi
if [ ! -f "$SLAM_BINDIR/$SLAM_BIN" ]; then							# Check SLAM bin if exists
	echo "SLAM bin does not exist in the specified path !"
	exit 1
fi
if [ -z $NCORE ]; then										# Check number of cores for calculations
	echo "Number of cores for calculations are not set, falling back to serial mode ..."
	NCORE="1"
fi
if [ !  -f sp_cluster_species.txt ]; then							# Check if potential file (sp_cluster_species.txt)
	echo "Potential information for running slam does not exist !"
	echo "Check if 'sp_cluster_species.txt' file is in the working directory"
	#exit 1
fi
if [ ! -d sp_cluster_parameter_src ]; then							# Check if basis function info exists
	echo "Basis functions information for running slam does not exits !"
	echo "Chekc if the directory 'sp_cluster_parameter_src' is in the working directory"
	#exit 1
fi


###############################################################################################
# MAIN ... 
###############################################################################################

declare -a config_cla		# ARR SAVING CLASSIC ION INFO
declare -a config_sp		# ARR SAVING SP ION INFO

TITLE=$( sed -n 1p geo.txt )	# SAVING COMMENT LINE WRITTEN IN 'geo.txt'
ATOMS=$( sed -n 2p geo.txt )	# SAVING NUMBER OF SPECIES 
ARR=( $ATOMS )			#
CLA_N=${ARR[0]}			#
SP_N=${ARR[1]}			#
TOTAL_ATOM_N=$( echo "$CLA_N + $SP_N" | bc )

for ((i=0; i<$CLA_N; i++)); do				# READ CLASSIC ION INFO
	LINE=$( echo "2 + $i + 1" | bc )		#	
	READ_LINE=$( sed -n "$LINE"p $TAR )		#
	SPLITER=( $READ_LINE )				#
	for ((j=0; j<5; j++)); do			#
		config_cla[$i*5+$j]=${SPLITER[$j]}	#
	done						#
done							# SAVED IN 1D ARR 'config_cla' STRIDE OF '5' PER EACH ATOM

for ((i=0; i<$SP_N; i++)); do				# READ SP ION INFO
	LINE=$( echo "2 + $CLA_N + $i + 1" | bc )	#
	READ_LINE=$( sed -n "$LINE"p $TAR )		#
	SPLITER=( $READ_LINE )				#
	for ((j=0; j<4; j++)); do			#
		config_sp[$i*4+$j]=${SPLITER[$j]}	#
	done						#
done							# SAVED IN 1D ARR 'config_sp' STRIDE OF '4' PER EACH ATOM


# PRE DEFINED FUNCTIONS FOR STDOUT CLASSIC (MM) AND SP (QM) ION INFO
function write_cla()	{
	local CLA_N="$1"
	shift
	local arr=("$@")
	local k=""
	for ((k=0; k<$CLA_N; k++)); do
		printf "%2s%2s%12.6f%12.6f%12.6f\n" ${config_cla[$k*5+0]} ${config_cla[$k*5+1]} \
			${config_cla[$k*5+2]} ${config_cla[$k*5+3]} ${config_cla[$k*5+4]}
	done
}
function write_sp()	{
	local SP_N="$1"
	shift
	local arr=("$@")
	local k=""	

	for ((k=0; k<$SP_N; k++)); do
		printf "%2s%14.6f%12.6f%12.6f\n" ${config_sp[$k*4+0]} ${config_sp[$k*4+1]} \
			${config_sp[$k*4+2]} ${config_sp[$k*4+3]}
	done
}
function header()	{
	echo $TITLE   >> $3
	echo "$1  $2" >> $3
}
# PRE DEFINED FUNCTIONS FOR STDOUT CLASSIC (MM) AND SP (QM) ION INFO ... END

printf "\n READ CHECK DONE ...\n"
printf "\n INPUT CONFIGURATION FROM 'geo.txt' \n"
printf "\n"
printf "%s\n" "--------------------------------------------------------------------------------------"
printf "\n"
write_cla $CLA_N ${config_cla[@]}
write_sp  $SP_N  ${config_sp[@]}
printf "\n"
printf "%s\n" "--------------------------------------------------------------------------------------"
printf "\n"

###############################################################################################
# GENERATE INPUT FILES FOR FDM

# SHELL ARR INDEX ... DO NOT MAKE ANY CHANGES !!!
CX="2"
CY="3"
CZ="4"
SX="1"
SY="2"
SZ="3"
FORCE_STR="Geometric Derivatives ( eV / Angstrom )"
# SHELL ARR INDEX ... DO NOT MAKE ANY CHANGES !!! END

mv geo.txt geo.txt.org
# GEN INPUT + RUN SINGLE POINTS
for ((i=0; i<"$SP_N+$CLA_N"; i++)); do
	INDEX=$(echo "$i+1" | bc)
	echo " WORKING ON ATOM $INDEX ..."
	printf "\n"
	# FDM CLA
	if [ $i -lt $CLA_N ]; then

		# PRE-PROCESS FOR CALCULATING HESSIAN ...  
		ATOM_TYPE_CLA=${config_cla[0]}
		ATOM_TYPE_SP=${config_sp[0]}
		ATOM_N=$( echo "$CLA_N + $SP_N" | bc )
		READ_N=$( echo "$ATOM_N + 6" | bc )
		printf "%4.d%12.d\n" $CLA_N $SP_N > atom_"$INDEX"_x

		# +X DELTA ... 9 lines
		ORIGINAL_C=${config_cla[$i*5+$CX]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"+x )
		FDM_C=$( echo "$ORIGINAL_C+$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_cla[$i*5+$CX]=$FDM_C
		write_cla $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp  $CLA_N ${config_sp[@]}  >> $FILE_NAME
		config_cla[$i*5+$CX]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz 0 > $FILE_NAME.out			#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		grep -A $READ_N "$FORCE_STR" $FILE_NAME.out > fdm_tmp
		grep "$ATOM_TYPE_CLA" fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$3,$4,$5}' >> atom_"$INDEX"_x
		grep "$ATOM_TYPE_SP"  fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$2,$3,$4}' >> atom_"$INDEX"_x
		echo "" >> atom_"$INDEX"_x
		rm fdm_tmp

		# -X DELTA
		ORIGINAL_C=${config_cla[$i*5+$CX]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"-x )
		FDM_C=$( echo "$ORIGINAL_C-$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_cla[$i*5+$CX]=$FDM_C
		write_cla $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp  $CLA_N ${config_sp[@]}  >> $FILE_NAME
		config_cla[$i*5+$CX]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz 0 > $FILE_NAME.out			#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#

		# RECORD FORCES
		grep -A $READ_N "$FORCE_STR" $FILE_NAME.out > fdm_tmp
		grep "$ATOM_TYPE_CLA" fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$3,$4,$5}' >> atom_"$INDEX"_x
		grep "$ATOM_TYPE_SP"  fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$2,$3,$4}' >> atom_"$INDEX"_x
		rm fdm_tmp

		# ----

		# PRE-PROCESS FOR CALCULATING HESSIAN ... 
		ATOM_TYPE_CLA=${config_cla[0]}
		ATOM_TYPE_SP=${config_sp[0]}
		ATOM_N=$( echo "$CLA_N + $SP_N" | bc )
		READ_N=$( echo "$ATOM_N + 6" | bc )
		printf "%4.d%12.d\n" $CLA_N $SP_N > atom_"$INDEX"_y

		# +Y DELTA ... 9 lines
		ORIGINAL_C=${config_cla[$i*5+$CY]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"+y )
		FDM_C=$( echo "$ORIGINAL_C+$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_cla[$i*5+$CY]=$FDM_C
		write_cla $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp  $CLA_N ${config_sp[@]}  >> $FILE_NAME
		config_cla[$i*5+$CY]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz 0 > $FILE_NAME.out			#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		grep -A $READ_N "$FORCE_STR" $FILE_NAME.out > fdm_tmp
		grep "$ATOM_TYPE_CLA" fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$3,$4,$5}' >> atom_"$INDEX"_y
		grep "$ATOM_TYPE_SP"  fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$2,$3,$4}' >> atom_"$INDEX"_y
		echo "" >> atom_"$INDEX"_y
		rm fdm_tmp

		# -Y DELTA
		ORIGINAL_C=${config_cla[$i*5+$CY]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"-y )
		FDM_C=$( echo "$ORIGINAL_C-$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_cla[$i*5+$CY]=$FDM_C
		write_cla $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp  $CLA_N ${config_sp[@]}  >> $FILE_NAME
		config_cla[$i*5+$CY]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz 0 > $FILE_NAME.out			#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		grep -A $READ_N "$FORCE_STR" $FILE_NAME.out > fdm_tmp
		grep "$ATOM_TYPE_CLA" fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$3,$4,$5}' >> atom_"$INDEX"_y
		grep "$ATOM_TYPE_SP"  fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$2,$3,$4}' >> atom_"$INDEX"_y
		rm fdm_tmp

		# ----

		# PRE-PROCESS FOR CALCULATING HESSIAN ... 
		ATOM_TYPE_CLA=${config_cla[0]}
		ATOM_TYPE_SP=${config_sp[0]}
		ATOM_N=$( echo "$CLA_N + $SP_N" | bc )
		READ_N=$( echo "$ATOM_N + 6" | bc )
		printf "%4.d%12.d\n" $CLA_N $SP_N > atom_"$INDEX"_z

		# +Z DELTA ... 9 lines
		ORIGINAL_C=${config_cla[$i*5+$CZ]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"+z )
		FDM_C=$( echo "$ORIGINAL_C+$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_cla[$i*5+$CZ]=$FDM_C
		write_cla $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp  $CLA_N ${config_sp[@]}  >> $FILE_NAME
		config_cla[$i*5+$CZ]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz 0 > $FILE_NAME.out			#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		grep -A $READ_N "$FORCE_STR" $FILE_NAME.out > fdm_tmp
		grep "$ATOM_TYPE_CLA" fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$3,$4,$5}' >> atom_"$INDEX"_z
		grep "$ATOM_TYPE_SP"  fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$2,$3,$4}' >> atom_"$INDEX"_z
		echo "" >> atom_"$INDEX"_z
		rm fdm_tmp

		# -Z DELTA
		ORIGINAL_C=${config_cla[$i*5+$CZ]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"-z )
		FDM_C=$( echo "$ORIGINAL_C-$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_cla[$i*5+$CZ]=$FDM_C
		write_cla $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp  $CLA_N ${config_sp[@]}  >> $FILE_NAME
		config_cla[$i*5+$CZ]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz 0 > $FILE_NAME.out			#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		grep -A $READ_N "$FORCE_STR" $FILE_NAME.out > fdm_tmp
		grep "$ATOM_TYPE_CLA" fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$3,$4,$5}' >> atom_"$INDEX"_z
		grep "$ATOM_TYPE_SP"  fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$2,$3,$4}' >> atom_"$INDEX"_z
		rm fdm_tmp

	else
	# FDM SP #####################################################################################################
		OFFSET=$( echo "$i-$CLA_N" | bc )

		# PRE-PROCESS FOR CALCULATING HESSIAN ...  
		ATOM_TYPE_CLA=${config_cla[0]}
		ATOM_TYPE_SP=${config_sp[0]}
		ATOM_N=$( echo "$CLA_N + $SP_N" | bc )
		READ_N=$( echo "$ATOM_N + 6" | bc )
		printf "%4.d%12.d\n" $CLA_N $SP_N > atom_"$INDEX"_x

		# +X DELTA ... 9 lines
		ORIGINAL_C=${config_sp[$OFFSET*4+$SX]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"+x )
		FDM_C=$( echo "$ORIGINAL_C+$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_sp[$OFFSET*4+$SX]=$FDM_C
		write_cla $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp  $CLA_N ${config_sp[@]}  >> $FILE_NAME
		config_sp[$OFFSET*4+$SX]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz 0 > $FILE_NAME.out			#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		grep -A $READ_N "$FORCE_STR" $FILE_NAME.out > fdm_tmp
		grep "$ATOM_TYPE_CLA" fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$3,$4,$5}' >> atom_"$INDEX"_x
		grep "$ATOM_TYPE_SP"  fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$2,$3,$4}' >> atom_"$INDEX"_x
		echo "" >> atom_"$INDEX"_x
		rm fdm_tmp

		# -X DELTA
		ORIGINAL_C=${config_sp[$OFFSET*4+$SX]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"-x )
		FDM_C=$( echo "$ORIGINAL_C-$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_sp[$OFFSET*4+$SX]=$FDM_C
		write_cla $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp  $CLA_N ${config_sp[@]}  >> $FILE_NAME
		config_sp[$OFFSET*4+$SX]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz 0 > $FILE_NAME.out			#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		grep -A $READ_N "$FORCE_STR" $FILE_NAME.out > fdm_tmp
		grep "$ATOM_TYPE_CLA" fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$3,$4,$5}' >> atom_"$INDEX"_x
		grep "$ATOM_TYPE_SP"  fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$2,$3,$4}' >> atom_"$INDEX"_x
		rm fdm_tmp

		# ----

		# PRE-PROCESS FOR CALCULATING HESSIAN ...  
		ATOM_TYPE_CLA=${config_cla[0]}
		ATOM_TYPE_SP=${config_sp[0]}
		ATOM_N=$( echo "$CLA_N + $SP_N" | bc )
		READ_N=$( echo "$ATOM_N + 6" | bc )
		printf "%4.d%12.d\n" $CLA_N $SP_N > atom_"$INDEX"_y

		# +Y DELTA ... 9 lines
		ORIGINAL_C=${config_sp[$OFFSET*4+$SY]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"+y )
		FDM_C=$( echo "$ORIGINAL_C+$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_sp[$OFFSET*4+$SY]=$FDM_C
		write_cla $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp  $CLA_N ${config_sp[@]}  >> $FILE_NAME
		config_sp[$OFFSET*4+$SY]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz 0 > $FILE_NAME.out			#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		grep -A $READ_N "$FORCE_STR" $FILE_NAME.out > fdm_tmp
		grep "$ATOM_TYPE_CLA" fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$3,$4,$5}' >> atom_"$INDEX"_y
		grep "$ATOM_TYPE_SP"  fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$2,$3,$4}' >> atom_"$INDEX"_y
		echo "" >> atom_"$INDEX"_y
		rm fdm_tmp

		# -Y DELTA
		ORIGINAL_C=${config_sp[$OFFSET*4+$SY]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"-y )
		FDM_C=$( echo "$ORIGINAL_C-$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_sp[$OFFSET*4+$SY]=$FDM_C
		write_cla $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp  $CLA_N ${config_sp[@]}  >> $FILE_NAME
		config_sp[$OFFSET*4+$SY]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz 0 > $FILE_NAME.out			#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		grep -A $READ_N "$FORCE_STR" $FILE_NAME.out > fdm_tmp
		grep "$ATOM_TYPE_CLA" fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$3,$4,$5}' >> atom_"$INDEX"_y
		grep "$ATOM_TYPE_SP"  fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$2,$3,$4}' >> atom_"$INDEX"_y
		rm fdm_tmp

		# ----

		# PRE-PROCESS FOR CALCULATING HESSIAN ...  
		ATOM_TYPE_CLA=${config_cla[0]}
		ATOM_TYPE_SP=${config_sp[0]}
		ATOM_N=$( echo "$CLA_N + $SP_N" | bc )
		READ_N=$( echo "$ATOM_N + 6" | bc )
		printf "%4.d%12.d\n" $CLA_N $SP_N > atom_"$INDEX"_z

		# +Z DELTA ... 9 lines
		ORIGINAL_C=${config_sp[$OFFSET*4+$SZ]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"+z )
		FDM_C=$( echo "$ORIGINAL_C+$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_sp[$OFFSET*4+$SZ]=$FDM_C
		write_cla $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp  $CLA_N ${config_sp[@]}  >> $FILE_NAME
		config_sp[$OFFSET*4+$SZ]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz 0 > $FILE_NAME.out			#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		grep -A $READ_N "$FORCE_STR" $FILE_NAME.out > fdm_tmp
		grep "$ATOM_TYPE_CLA" fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$3,$4,$5}' >> atom_"$INDEX"_z
		grep "$ATOM_TYPE_SP"  fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$2,$3,$4}' >> atom_"$INDEX"_z
		echo "" >> atom_"$INDEX"_z
		rm fdm_tmp

		# -Z DELTA
		ORIGINAL_C=${config_sp[$OFFSET*4+$SZ]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"-z )
		FDM_C=$( echo "$ORIGINAL_C-$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_sp[$OFFSET*4+$SZ]=$FDM_C
		write_cla $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp  $CLA_N ${config_sp[@]}  >> $FILE_NAME
		config_sp[$OFFSET*4+$SZ]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz 0 > $FILE_NAME.out			#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		grep -A $READ_N "$FORCE_STR" $FILE_NAME.out > fdm_tmp
		grep "$ATOM_TYPE_CLA" fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$3,$4,$5}' >> atom_"$INDEX"_z
		grep "$ATOM_TYPE_SP"  fdm_tmp | awk '{printf " %3s%16s%16s%16s\n", $1,$2,$3,$4}' >> atom_"$INDEX"_z
		rm fdm_tmp
	fi

done
mv geo.txt.org geo.txt

printf "%s\n" " EXECUTING NORMAL MODE ANALYSIS ..."
printf "\n"

python $DIGPY $TOTAL_ATOM_N $DELTA			# CALL FREQ CALCULATOR ... USING WILSON'S GF METHOD FOR NORMAL MODE CALC
bash   $CLEAN
printf "%s\n" "######################################################################################"
