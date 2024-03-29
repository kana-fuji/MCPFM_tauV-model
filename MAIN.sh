#!/bin/bash
#---run
#############################################################################
#

#
#
#---compile
#############################################################################
#

#
#
#---func
#############################################################################

fmkdir () {
    if [ ! -d $1 ]; then
	mkdir $1
    fi
}


#
#---main
#############################################################################

fmkdir "DATA"

Dir="DATA/"$1
fmkdir $Dir
rm -fr $Dir/*

fmkdir $Dir"/INPUT"
fmkdir $Dir"/Code"
fmkdir $Dir"/Debug"
cp * ${Dir}/Code

date> ${Dir}/time
start_time=`date +%s`

./run $1 $2 $3 $4 $5 $6 $7 $8 $9 > ${Dir}/out
#./mcpf_2d_us $1 < IN/$1.param > ${Dir}/out

date>> ${Dir}/time
end_time=`date +%s`

RT=`expr ${end_time} - ${start_time}`
H=`expr ${RT} / 3600`
RT=`expr ${RT} % 3600`
M=`expr ${RT} / 60`
S=`expr ${RT} % 60`
echo "${H}:${M}:${S}" >> ${Dir}/time


