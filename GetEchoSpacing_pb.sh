#!/bin/bash
PURPOSE="Compute B0 Unwarping factors from a .PAR file" ;



#1.0 6/29/10	-- CJ
#1.1 7/15/10 	-- Add input error detection routines


USAGE() {
    cat <<EOF

	 ===+++===+++     Program:  ${PGM}
	 Use this routine to:
	 	${PURPOSE}
	 Usage: 
	  GetEchoSpacing <MyFile.PAR> <SENSE Factor>
	 
	 where :
	   <MyFile.PAR> is an fMRI PAR file
	   <SENSE Factor> is your SENSE Acceleration Factor from the exam card
	 	 
EOF

	}


#
#
# test number and value of arguments

#clear
EXITCODE=0
if [ $# -lt 2 ] ; then
	# called without input file?
	EXITCODE=1
elif [ ! -e ${1} ] ; then
	# input file does not exist
	echo 
	echo
	echo The file: ${1}
	echo Does not exist in the folder `pwd`
	echo
	EXITCODE=1
else
	# Is this a PAR file?

	PARFILETST=`echo ${1} | awk '{ print index($0,".PAR")}'` ;
#echo $PARFILETST ;

	if [   ${PARFILETST} -eq 0   ]; then
		echo ${1} is not a .PAR file
		EXITCODE=1
	fi
fi


if [ ${EXITCODE} -eq 1 ] ; then
	USAGE ;
	echo
	echo PAR files in this folder are listed below:
	echo
	ls *.PAR
	exit
fi


THEPARFILE=${1}
ACCELERATION=${2}
DIR=`dirname $THEPARFILE`
BASE=`basename $THEPARFILE .PAR`

# Collect data from the Input file

#Dynamic scan      <0=no 1=yes> ?   :   1
	DYN=`awk '/Dynamic scan/ {print $8}' ${THEPARFILE}`
#	echo Dynamic scan: $DYN

#.    Diffusion         <0=no 1=yes> ?   :   0
	DIFF=`awk '/.    Diffusion     / {print $7}' ${THEPARFILE}`
#	echo Diffusion: $DIFF
#.    Number of label types   <0=no ASL> :   0
	LABELTYPES=`awk '/.    Number of label types/ {print int($9)}' ${THEPARFILE}`
#	echo Number of label types: $LABELTYPES
#.    Water Fat shift [pixels]           :   10.881
	WFS=`awk '/Water Fat shift/ {print $7}' ${THEPARFILE}`
#	echo Water Fat shift: ${WFS}
#.    EPI factor        <0,1=no EPI>     :   63
	EPIFACTOR=`awk '/EPI factor/ {print $7}'  ${THEPARFILE}`
#	echo EPI factor: ${EPIFACTOR}
# First slice 
	ECHO=`awk '/  1   1    1  / {print $31}' ${THEPARFILE}`
#	echo Echo: ${ECHO}
# test if this is an fMRI PAR file
# fMRI PAR file iff
#	EPIFACTOR > 1
#	DYN = 1
#	DIFF = 0
#	LABELTYPES > 0

EXITCODE=0;
if [ ${EPIFACTOR} \< 2 ] ; then
echo Not an EPI protocol
EXITCODE=1
fi
if [ ${DYN} = 0 ] ; then
echo Not a Dynamic scan
EXITCODE=1
fi
if [ ${DIFF} = 1 ] ; then
echo Diffusion Protocol
EXITCODE=1
fi
# can't get the following to work so just skip it
if [ ${LABELTYPES} !=  0 ] ; then
echo LABELTYPES indicate an ASL protocol
EXITCODE=1
fi

#echo ExitCode: ${EXITCODE}
#echo
if [ ${EXITCODE} = 1 ] ; then
echo Input file: $THEPARFILE is not an fMRI PAR file
echo
echo
exit
fi

#echo Input file: $THEPARFILE
# Compute
#echo
python << EOF
WFS = ${WFS}
EPIFACTOR = ${EPIFACTOR}
ACCELERATION = ${ACCELERATION}
DIR = '${DIR}'
BASE = '${BASE}'
FACTOR = (1000*WFS)/(434.215*(EPIFACTOR + 1))/(ACCELERATION*1000)
#print 'EchoSpacing = %.8f seconds' % FACTOR
print '%.8f' % FACTOR
#text_file = open((DIR + '/' + BASE + '_echo_spacing.txt'),"w")
#text_file.write('%.8f' % FACTOR)
#text_file.close()

EOF
#echo $FACTOR

#echo
#echo

