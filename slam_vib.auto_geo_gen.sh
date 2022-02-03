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
printf "%s\n" " compatibilty : above SLAM V 2.2.1"
printf "%s\n" " last edited  : 31 - 08 -2020"
printf "%s\n" "--------------------------------------------------------------------------------------"

# Default env - variables.
#OPTI_MX=5000
TAR="geo.txt"								# Takes standard slam input geometry file 'geo.txt'
DELTA=$1								# Takes "dx" value, if it is not specified. then  ... use defualt value below

printf "\n"
if [ -z $DELTA ]; then							# falling back to the default value '0.0025'
	DELTA="0.0025"							# default delta value ... 
	printf "%2s %12.6f\n" " Default 'dx' value is used :" $DELTA	#
else									#
	printf "%2s %12.6f\n" " Requested 'dx' value is used :" $DELTA	#
fi									#

# SOME Environmental variables // user must set those depending on your setup

# SLAM bin path & excutable ... user customisation available
#SLAM_BINDIR="/home/uccawkj/src_tool/bin"				# SLAM BIN DIR
#SLAM_BIN="slam_v_2.2_opti.mpi.x.bipb"					# SLAM BIN PATH
SLAM_BINDIR="/home/uccawkj/bin"						# SLAM BIN DIR
SLAM_BIN="slam.240122.mpi.x"					# SLAM BIN PATH
PREFIX="mpirun -np"							# Parallel process arguments
NCORE="2"								# Number of cores requested for the calculation (single point)
EXE="$PREFIX $NCORE $SLAM_BINDIR/$SLAM_BIN"
DIGPY="/home/uccawkj/bin/slam_vib_package_v_2.2/slam_vib.diagonaliser.py"
CLEAN="/home/uccawkj/bin/slam_vib_package_v_2.2/clean.sh"
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
	echo "Basis function information for running slam does not exits !"
	echo "Check if the directory 'sp_cluster_parameter_src' is in the working directory"
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
CLA_N=${ARR[0]}			# GET CLASSICAL SPECIES #
SP_N=${ARR[1]}			# GET SP SPECIES #
#TOTAL_ATOM_N=$( echo "$CLA_N + $SP_N" | bc )

for ((i=0; i<$CLA_N; i++)); do				# READ CLASSIC ION INFO
	LINE=$( echo "2 + $i + 1" | bc )		#	
	READ_LINE=$( sed -n "$LINE"p $TAR )		#
	SPLITER=( $READ_LINE )				#
	for ((j=0; j<5; j++)); do			#
		config_cla[$i*5+$j]=${SPLITER[$j]}	# STRIDE 5 ... name kind(c/s) x y z
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
function write_cla_parser()	{
	local CLA_N="$1"
	shift
	local arr=("$@")
	local k=""
	for ((k=0; k<$CLA_N; k++)); do

		if [ ${config_cla[$k*5+1]} == 'c' ]; then
			printf "%2s%2s%12.6f%12.6f%12.6f%12s\n" ${config_cla[$k*5+0]} ${config_cla[$k*5+1]} \
				${config_cla[$k*5+2]} ${config_cla[$k*5+3]} ${config_cla[$k*5+4]} "_fix_"
		else
			printf "%2s%2s%12.6f%12.6f%12.6f\n" ${config_cla[$k*5+0]} ${config_cla[$k*5+1]} \
				${config_cla[$k*5+2]} ${config_cla[$k*5+3]} ${config_cla[$k*5+4]}
		fi
	done
}
function write_sp_parser()	{
	local SP_N="$1"
	shift
	local arr=("$@")
	local k=""	

	for ((k=0; k<$SP_N; k++)); do
		printf "%2s%14.6f%12.6f%12.6f%12s\n" ${config_sp[$k*4+0]} ${config_sp[$k*4+1]} \
			${config_sp[$k*4+2]} ${config_sp[$k*4+3]} "_fix_"
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


# CLASSIC SHELL MODEL CHEKC ... IF THERE IS NO SHELL THEN DO PURE SINGLE CALC !!! >>> OPTI_MX = 0

if_shell_used=0							# IF SHELL CLASSIC ION USED CHEKCER FLAG
for ((i=0; i<$CLA_N; i++)); do					#
	if [ ${config_cla[$i*5+1]} == 's' ]; then		#
		if_shell_used=1					# IF SHELL USED ... SET FLAG TRUE
		break		
	fi
done
if [ $if_shell_used == 1 ]; then				# IF SHELL USED ... SET OPTI SCF CYCLE WITH DEFAULRT VALUE 5000
#DEBUG
	OPTI_MX=5000						#
else								#
	OPTI_MX=0						# ELSE DO SINGLE ONLY
fi								#


mv geo.txt geo.txt.org
# GEN INPUT + RUN SINGLE POINTS
INDEX=0

for ((i=0; i<"$SP_N+$CLA_N"; i++)); do

	# FDM CLA
	#if [ $i -lt $CLA_N ] && [ ${config_cla[$i*5+1]} == 'c' ]; then		# CHECK IF THE ATOM IS CORE
	if [ $i -lt $CLA_N ]; then

		if [ ${config_cla[$i*5+1]} == 'c' ]; then

		INDEX=$(echo "$INDEX + 1" | bc )
		echo " WORKING ON ATOM $INDEX ..."
		printf "\n"

		# PRE-PROCESS FOR CALCULATING HESSIAN ...  
		ATOM_N=$( echo "$CLA_N + $SP_N" | bc )
		READ_N=$( echo "$ATOM_N + 4" | bc )

		# +X DELTA ... 9 lines
		ORIGINAL_C=${config_cla[$i*5+$CX]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"+x )
		FDM_C=$( echo "$ORIGINAL_C+$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_cla[$i*5+$CX]=$FDM_C
		write_cla_parser $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp_parser  $SP_N  ${config_sp[@]}  >> $FILE_NAME
		config_cla[$i*5+$CX]=$ORIGINAL_C			# PUT BACK THE ORIGINAL INPUT COORDINATE
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz $OPTI_MX > $FILE_NAME.out		#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
	
		
		# RECORD FORCES
		line=$(grep -n "Geometric Derivatives ( eV / Angstrom )" "$FILE_NAME".out | awk '{print $1}' | tail -1)
		line=${line:0:-1}
		sta_line=$( echo "$line + 5" | bc )
		end_line=$( echo "$line + $READ_N" | bc )
		sed -n "$sta_line","$end_line"p $FILE_NAME.out  > fdm_tmp

		CLA_CORE_N=$(grep " c " fdm_tmp | wc -l | awk '{print $1}')

		printf "%4.d%12.d\n" $CLA_CORE_N $SP_N > atom_"$INDEX"_x
		for (( k=1; k<=$CLA_N + $SP_N; k++ )); do
			rl=$(sed -n "$k"p fdm_tmp)
			spl=( $rl )
			if [ $k -le $CLA_N ]; then
				if [ ${spl[1]} == 'c' ]; then
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[2]} ${spl[3]} ${spl[4]} >> atom_"$INDEX"_x
				fi
			else
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[1]} ${spl[2]} ${spl[3]} >> atom_"$INDEX"_x
			fi
		done

		echo "" >> atom_"$INDEX"_x
		rm fdm_tmp


		# -X DELTA ... 9 lines
		ORIGINAL_C=${config_cla[$i*5+$CX]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"-x )
		FDM_C=$( echo "$ORIGINAL_C-$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_cla[$i*5+$CX]=$FDM_C
		write_cla_parser $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp_parser  $SP_N  ${config_sp[@]}  >> $FILE_NAME
		config_cla[$i*5+$CX]=$ORIGINAL_C			# PUT BACK THE ORIGINAL INPUT COORDINATE
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz $OPTI_MX > $FILE_NAME.out		#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
	
		
		# RECORD FORCES
		line=$(grep -n "Geometric Derivatives ( eV / Angstrom )" "$FILE_NAME".out | awk '{print $1}' | tail -1)
		line=${line:0:-1}
		sta_line=$( echo "$line + 5" | bc )
		end_line=$( echo "$line + $READ_N" | bc )
		sed -n "$sta_line","$end_line"p $FILE_NAME.out  > fdm_tmp

		CLA_CORE_N=$(grep " c " fdm_tmp | wc -l | awk '{print $1}')

		#printf "%4.d%12.d\n" $CLA_CORE_N $SP_N > atom_"$INDEX"_x
		for (( k=1; k<=$CLA_N + $SP_N; k++ )); do
			rl=$(sed -n "$k"p fdm_tmp)
			spl=( $rl )
			if [ $k -le $CLA_N ]; then
				if [ ${spl[1]} == 'c' ]; then
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[2]} ${spl[3]} ${spl[4]} >> atom_"$INDEX"_x
				fi
			else
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[1]} ${spl[2]} ${spl[3]} >> atom_"$INDEX"_x
			fi
		done

		#echo "" >> atom_"$INDEX"_x
		rm fdm_tmp

		# ----

		# PRE-PROCESS FOR CALCULATING HESSIAN ...  
		ATOM_N=$( echo "$CLA_N + $SP_N" | bc )
		READ_N=$( echo "$ATOM_N + 4" | bc )

		# +Y DELTA ... 9 lines
		ORIGINAL_C=${config_cla[$i*5+$CY]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"+y )
		FDM_C=$( echo "$ORIGINAL_C+$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_cla[$i*5+$CY]=$FDM_C
		write_cla_parser $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp_parser  $SP_N  ${config_sp[@]}  >> $FILE_NAME
		config_cla[$i*5+$CY]=$ORIGINAL_C			# PUT BACK THE ORIGINAL INPUT COORDINATE
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz $OPTI_MX > $FILE_NAME.out		#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
	
		
		# RECORD FORCES
		line=$(grep -n "Geometric Derivatives ( eV / Angstrom )" "$FILE_NAME".out | awk '{print $1}' | tail -1)
		line=${line:0:-1}
		sta_line=$( echo "$line + 5" | bc )
		end_line=$( echo "$line + $READ_N" | bc )
		sed -n "$sta_line","$end_line"p $FILE_NAME.out  > fdm_tmp

		CLA_CORE_N=$(grep " c " fdm_tmp | wc -l | awk '{print $1}')

		printf "%4.d%12.d\n" $CLA_CORE_N $SP_N > atom_"$INDEX"_y
		for (( k=1; k<=$CLA_N + $SP_N; k++ )); do
			rl=$(sed -n "$k"p fdm_tmp)
			spl=( $rl )
			if [ $k -le $CLA_N ]; then
				if [ ${spl[1]} == 'c' ]; then
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[2]} ${spl[3]} ${spl[4]} >> atom_"$INDEX"_y
				fi
			else
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[1]} ${spl[2]} ${spl[3]} >> atom_"$INDEX"_y
			fi
		done

		echo "" >> atom_"$INDEX"_y
		rm fdm_tmp


		# -Y DELTA ... 9 lines
		ORIGINAL_C=${config_cla[$i*5+$CY]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"-y )
		FDM_C=$( echo "$ORIGINAL_C-$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_cla[$i*5+$CY]=$FDM_C
		write_cla_parser $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp_parser  $SP_N  ${config_sp[@]}  >> $FILE_NAME
		config_cla[$i*5+$CY]=$ORIGINAL_C			# PUT BACK THE ORIGINAL INPUT COORDINATE
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz $OPTI_MX > $FILE_NAME.out		#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
	
		
		# RECORD FORCES
		line=$(grep -n "Geometric Derivatives ( eV / Angstrom )" "$FILE_NAME".out | awk '{print $1}' | tail -1)
		line=${line:0:-1}
		sta_line=$( echo "$line + 5" | bc )
		end_line=$( echo "$line + $READ_N" | bc )
		sed -n "$sta_line","$end_line"p $FILE_NAME.out  > fdm_tmp

		CLA_CORE_N=$(grep " c " fdm_tmp | wc -l | awk '{print $1}')

		#printf "%4.d%12.d\n" $CLA_CORE_N $SP_N > atom_"$INDEX"_x
		for (( k=1; k<=$CLA_N + $SP_N; k++ )); do
			rl=$(sed -n "$k"p fdm_tmp)
			spl=( $rl )
			if [ $k -le $CLA_N ]; then
				if [ ${spl[1]} == 'c' ]; then
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[2]} ${spl[3]} ${spl[4]} >> atom_"$INDEX"_y
				fi
			else
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[1]} ${spl[2]} ${spl[3]} >> atom_"$INDEX"_y
			fi
		done

		#echo "" >> atom_"$INDEX"_y
		rm fdm_tmp

		# ----

		# PRE-PROCESS FOR CALCULATING HESSIAN ...  
		ATOM_N=$( echo "$CLA_N + $SP_N" | bc )
		READ_N=$( echo "$ATOM_N + 4" | bc )

		# +Z DELTA ... 9 lines
		ORIGINAL_C=${config_cla[$i*5+$CZ]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"+z )
		FDM_C=$( echo "$ORIGINAL_C+$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_cla[$i*5+$CZ]=$FDM_C
		write_cla_parser $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp_parser  $SP_N  ${config_sp[@]}  >> $FILE_NAME
		config_cla[$i*5+$CZ]=$ORIGINAL_C			# PUT BACK THE ORIGINAL INPUT COORDINATE
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz $OPTI_MX > $FILE_NAME.out		#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
	
		
		# RECORD FORCES
		line=$(grep -n "Geometric Derivatives ( eV / Angstrom )" "$FILE_NAME".out | awk '{print $1}' | tail -1)
		line=${line:0:-1}
		sta_line=$( echo "$line + 5" | bc )
		end_line=$( echo "$line + $READ_N" | bc )
		sed -n "$sta_line","$end_line"p $FILE_NAME.out  > fdm_tmp

		CLA_CORE_N=$(grep " c " fdm_tmp | wc -l | awk '{print $1}')

		printf "%4.d%12.d\n" $CLA_CORE_N $SP_N > atom_"$INDEX"_z
		for (( k=1; k<=$CLA_N + $SP_N; k++ )); do
			rl=$(sed -n "$k"p fdm_tmp)
			spl=( $rl )
			if [ $k -le $CLA_N ]; then
				if [ ${spl[1]} == 'c' ]; then
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[2]} ${spl[3]} ${spl[4]} >> atom_"$INDEX"_z
				fi
			else
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[1]} ${spl[2]} ${spl[3]} >> atom_"$INDEX"_z
			fi
		done

		echo "" >> atom_"$INDEX"_z
		rm fdm_tmp


		# -Z DELTA ... 9 lines
		ORIGINAL_C=${config_cla[$i*5+$CZ]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"-z )
		FDM_C=$( echo "$ORIGINAL_C-$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_cla[$i*5+$CZ]=$FDM_C
		write_cla_parser $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp_parser  $SP_N  ${config_sp[@]}  >> $FILE_NAME
		config_cla[$i*5+$CZ]=$ORIGINAL_C			# PUT BACK THE ORIGINAL INPUT COORDINATE
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz $OPTI_MX > $FILE_NAME.out		#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
	
		
		# RECORD FORCES
		line=$(grep -n "Geometric Derivatives ( eV / Angstrom )" "$FILE_NAME".out | awk '{print $1}' | tail -1)
		line=${line:0:-1}
		sta_line=$( echo "$line + 5" | bc )
		end_line=$( echo "$line + $READ_N" | bc )
		sed -n "$sta_line","$end_line"p $FILE_NAME.out  > fdm_tmp

		CLA_CORE_N=$(grep " c " fdm_tmp | wc -l | awk '{print $1}')

		#printf "%4.d%12.d\n" $CLA_CORE_N $SP_N > atom_"$INDEX"_z
		for (( k=1; k<=$CLA_N + $SP_N; k++ )); do
			rl=$(sed -n "$k"p fdm_tmp)
			spl=( $rl )
			if [ $k -le $CLA_N ]; then
				if [ ${spl[1]} == 'c' ]; then
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[2]} ${spl[3]} ${spl[4]} >> atom_"$INDEX"_z
				fi
			else
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[1]} ${spl[2]} ${spl[3]} >> atom_"$INDEX"_z
			fi
		done

		#echo "" >> atom_"$INDEX"_z
		rm fdm_tmp
	
		fi
		# ----

	else
	# FDM SP #####################################################################################################
		OFFSET=$( echo "$i-$CLA_N" | bc )

		# INDEX UPDATE
		INDEX=$(echo "$INDEX + 1" | bc )
		echo " WORKING ON ATOM $INDEX ..."
		printf "\n"

		# PRE-PROCESS FOR CALCULATING HESSIAN ...  
		ATOM_N=$( echo "$CLA_N + $SP_N" | bc )
		READ_N=$( echo "$ATOM_N + 4" | bc )

		# +X DELTA ... 9 lines
		ORIGINAL_C=${config_sp[$OFFSET*4+$SX]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"+x )
		FDM_C=$( echo "$ORIGINAL_C+$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_sp[$OFFSET*4+$SX]=$FDM_C
		write_cla_parser $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp_parser  $SP_N ${config_sp[@]}  >> $FILE_NAME
		config_sp[$OFFSET*4+$SX]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz $OPTI_MX > $FILE_NAME.out		#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		line=$(grep -n "Geometric Derivatives ( eV / Angstrom )" "$FILE_NAME".out | awk '{print $1}' | tail -1)
		line=${line:0:-1}
		sta_line=$( echo "$line + 5" | bc )
		end_line=$( echo "$line + $READ_N" | bc )
		sed -n "$sta_line","$end_line"p $FILE_NAME.out  > fdm_tmp

		CLA_CORE_N=$(grep " c " fdm_tmp | wc -l | awk '{print $1}')

		printf "%4.d%12.d\n" $CLA_CORE_N $SP_N > atom_"$INDEX"_x
		for (( k=1; k<=$CLA_N + $SP_N; k++ )); do
			rl=$(sed -n "$k"p fdm_tmp)
			spl=( $rl )
			if [ $k -le $CLA_N ]; then
				if [ ${spl[1]} == 'c' ]; then
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[2]} ${spl[3]} ${spl[4]} >> atom_"$INDEX"_x
				fi
			else
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[1]} ${spl[2]} ${spl[3]} >> atom_"$INDEX"_x
			fi
		done

		echo "" >> atom_"$INDEX"_x
		rm fdm_tmp

		# PRE-PROCESS FOR CALCULATING HESSIAN ...  
		ATOM_N=$( echo "$CLA_N + $SP_N" | bc )
		READ_N=$( echo "$ATOM_N + 4" | bc )

		# -X DELTA ... 9 lines
		ORIGINAL_C=${config_sp[$OFFSET*4+$SX]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"-x )
		FDM_C=$( echo "$ORIGINAL_C-$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_sp[$OFFSET*4+$SX]=$FDM_C
		write_cla_parser $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp_parser  $SP_N ${config_sp[@]}  >> $FILE_NAME
		config_sp[$OFFSET*4+$SX]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz $OPTI_MX > $FILE_NAME.out		#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		line=$(grep -n "Geometric Derivatives ( eV / Angstrom )" "$FILE_NAME".out | awk '{print $1}' | tail -1)
		line=${line:0:-1}
		sta_line=$( echo "$line + 5" | bc )
		end_line=$( echo "$line + $READ_N" | bc )
		sed -n "$sta_line","$end_line"p $FILE_NAME.out  > fdm_tmp

		CLA_CORE_N=$(grep " c " fdm_tmp | wc -l | awk '{print $1}')

		#printf "%4.d%12.d\n" $CLA_CORE_N $SP_N > atom_"$INDEX"_x
		for (( k=1; k<=$CLA_N + $SP_N; k++ )); do
			rl=$(sed -n "$k"p fdm_tmp)
			spl=( $rl )
			if [ $k -le $CLA_N ]; then
				if [ ${spl[1]} == 'c' ]; then
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[2]} ${spl[3]} ${spl[4]} >> atom_"$INDEX"_x
				fi
			else
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[1]} ${spl[2]} ${spl[3]} >> atom_"$INDEX"_x
			fi
		done

		rm fdm_tmp

		# ------------

		# PRE-PROCESS FOR CALCULATING HESSIAN ...  
		ATOM_N=$( echo "$CLA_N + $SP_N" | bc )
		READ_N=$( echo "$ATOM_N + 4" | bc )

		# +Y DELTA ... 9 lines
		ORIGINAL_C=${config_sp[$OFFSET*4+$SY]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"+y )
		FDM_C=$( echo "$ORIGINAL_C+$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_sp[$OFFSET*4+$SY]=$FDM_C
		write_cla_parser $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp_parser  $SP_N ${config_sp[@]}  >> $FILE_NAME
		config_sp[$OFFSET*4+$SY]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz $OPTI_MX > $FILE_NAME.out		#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		line=$(grep -n "Geometric Derivatives ( eV / Angstrom )" "$FILE_NAME".out | awk '{print $1}' | tail -1)
		line=${line:0:-1}
		sta_line=$( echo "$line + 5" | bc )
		end_line=$( echo "$line + $READ_N" | bc )
		sed -n "$sta_line","$end_line"p $FILE_NAME.out  > fdm_tmp

		CLA_CORE_N=$(grep " c " fdm_tmp | wc -l | awk '{print $1}')

		printf "%4.d%12.d\n" $CLA_CORE_N $SP_N > atom_"$INDEX"_y
		for (( k=1; k<=$CLA_N + $SP_N; k++ )); do
			rl=$(sed -n "$k"p fdm_tmp)
			spl=( $rl )
			if [ $k -le $CLA_N ]; then
				if [ ${spl[1]} == 'c' ]; then
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[2]} ${spl[3]} ${spl[4]} >> atom_"$INDEX"_y
				fi
			else
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[1]} ${spl[2]} ${spl[3]} >> atom_"$INDEX"_y
			fi
		done

		echo "" >> atom_"$INDEX"_y
		rm fdm_tmp

		# PRE-PROCESS FOR CALCULATING HESSIAN ...  
		ATOM_N=$( echo "$CLA_N + $SP_N" | bc )
		READ_N=$( echo "$ATOM_N + 4" | bc )

		# -Y DELTA ... 9 lines
		ORIGINAL_C=${config_sp[$OFFSET*4+$SY]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"-y )
		FDM_C=$( echo "$ORIGINAL_C-$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_sp[$OFFSET*4+$SY]=$FDM_C
		write_cla_parser $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp_parser  $SP_N ${config_sp[@]}  >> $FILE_NAME
		config_sp[$OFFSET*4+$SY]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz $OPTI_MX > $FILE_NAME.out		#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		line=$(grep -n "Geometric Derivatives ( eV / Angstrom )" "$FILE_NAME".out | awk '{print $1}' | tail -1)
		line=${line:0:-1}
		sta_line=$( echo "$line + 5" | bc )
		end_line=$( echo "$line + $READ_N" | bc )
		sed -n "$sta_line","$end_line"p $FILE_NAME.out  > fdm_tmp

		CLA_CORE_N=$(grep " c " fdm_tmp | wc -l | awk '{print $1}')

		#printf "%4.d%12.d\n" $CLA_CORE_N $SP_N > atom_"$INDEX"_y
		for (( k=1; k<=$CLA_N + $SP_N; k++ )); do
			rl=$(sed -n "$k"p fdm_tmp)
			spl=( $rl )
			if [ $k -le $CLA_N ]; then
				if [ ${spl[1]} == 'c' ]; then
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[2]} ${spl[3]} ${spl[4]} >> atom_"$INDEX"_y
				fi
			else
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[1]} ${spl[2]} ${spl[3]} >> atom_"$INDEX"_y
			fi
		done

		rm fdm_tmp

		# ------------

		# PRE-PROCESS FOR CALCULATING HESSIAN ...  
		ATOM_N=$( echo "$CLA_N + $SP_N" | bc )
		READ_N=$( echo "$ATOM_N + 4" | bc )

		# +Z DELTA ... 9 lines
		ORIGINAL_C=${config_sp[$OFFSET*4+$SZ]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"+z )
		FDM_C=$( echo "$ORIGINAL_C+$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_sp[$OFFSET*4+$SZ]=$FDM_C
		write_cla_parser $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp_parser  $SP_N ${config_sp[@]}  >> $FILE_NAME
		config_sp[$OFFSET*4+$SZ]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz $OPTI_MX > $FILE_NAME.out		#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		line=$(grep -n "Geometric Derivatives ( eV / Angstrom )" "$FILE_NAME".out | awk '{print $1}' | tail -1)
		line=${line:0:-1}
		sta_line=$( echo "$line + 5" | bc )
		end_line=$( echo "$line + $READ_N" | bc )
		sed -n "$sta_line","$end_line"p $FILE_NAME.out  > fdm_tmp

		CLA_CORE_N=$(grep " c " fdm_tmp | wc -l | awk '{print $1}')

		printf "%4.d%12.d\n" $CLA_CORE_N $SP_N > atom_"$INDEX"_z
		for (( k=1; k<=$CLA_N + $SP_N; k++ )); do
			rl=$(sed -n "$k"p fdm_tmp)
			spl=( $rl )
			if [ $k -le $CLA_N ]; then
				if [ ${spl[1]} == 'c' ]; then
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[2]} ${spl[3]} ${spl[4]} >> atom_"$INDEX"_z
				fi
			else
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[1]} ${spl[2]} ${spl[3]} >> atom_"$INDEX"_z
			fi
		done

		echo "" >> atom_"$INDEX"_z
		rm fdm_tmp

		# PRE-PROCESS FOR CALCULATING HESSIAN ...  
		ATOM_N=$( echo "$CLA_N + $SP_N" | bc )
		READ_N=$( echo "$ATOM_N + 4" | bc )

		# -Z DELTA ... 9 lines
		ORIGINAL_C=${config_sp[$OFFSET*4+$SZ]}
		FILE_NAME=$( echo geo.txt."$INDEX"_"$DELTA"-z )
		FDM_C=$( echo "$ORIGINAL_C-$DELTA" | bc -l )
		touch $FILE_NAME
		header $CLA_N $SP_N $FILE_NAME
		config_sp[$OFFSET*4+$SZ]=$FDM_C
		write_cla_parser $CLA_N ${config_cla[@]} >> $FILE_NAME
		write_sp_parser  $SP_N ${config_sp[@]}  >> $FILE_NAME
		config_sp[$OFFSET*4+$SZ]=$ORIGINAL_C
		mv $FILE_NAME geo.txt					# RUN SINLGE POINT CALC
		$EXE $FILE_NAME.xyz $OPTI_MX > $FILE_NAME.out		#
		rm geo.txt.next						#
		mv geo.txt $FILE_NAME					#
		
		# RECORD FORCES
		line=$(grep -n "Geometric Derivatives ( eV / Angstrom )" "$FILE_NAME".out | awk '{print $1}' | tail -1)
		line=${line:0:-1}
		sta_line=$( echo "$line + 5" | bc )
		end_line=$( echo "$line + $READ_N" | bc )
		sed -n "$sta_line","$end_line"p $FILE_NAME.out  > fdm_tmp

		CLA_CORE_N=$(grep " c " fdm_tmp | wc -l | awk '{print $1}')

		#printf "%4.d%12.d\n" $CLA_CORE_N $SP_N > atom_"$INDEX"_x
		for (( k=1; k<=$CLA_N + $SP_N; k++ )); do
			rl=$(sed -n "$k"p fdm_tmp)
			spl=( $rl )
			if [ $k -le $CLA_N ]; then
				if [ ${spl[1]} == 'c' ]; then
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[2]} ${spl[3]} ${spl[4]} >> atom_"$INDEX"_z
				fi
			else
					printf "%3s%16s%16s%16s\n" ${spl[0]} ${spl[1]} ${spl[2]} ${spl[3]} >> atom_"$INDEX"_z
			fi
		done

		rm fdm_tmp

		# ------------
	fi

done
mv geo.txt.org geo.txt

printf "%s\n" " EXECUTING NORMAL MODE ANALYSIS ..."
printf "\n"

TOTAL_ATOM_N=$( echo "$CLA_CORE_N + $SP_N" | bc )
python $DIGPY $TOTAL_ATOM_N $DELTA			# CALL FREQ CALCULATOR ... USING WILSON'S GF METHOD FOR NORMAL MODE CALC
bash   $CLEAN
printf "%s\n" "######################################################################################"
