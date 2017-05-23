#!/bin/bash
funcName=rest
physioName=${funcName}_physio
subjectsDir=/NAS_II/Projects/Udall/subjects
cwd=$(pwd)
for sessionDir in $(ls -d RC4???-?); do
	subject=$(echo $sessionDir|cut -d'-' -f 1)
	sessionNr=$(echo $sessionDir|cut -d'-' -f 2)
	cd $cwd/$sessionDir
	echo $(pwd)
	physioPath=$subjectsDir/$subject/session$sessionNr/physio/${physioName}.log
	if [ -e $physioPath ]; then
		#echo $physioPath
		#ln -sf $physioPath ${physioName}.log
		#sls_extractphysio ${physioName}.log 720 > ${physioName}_extracted.txt
		#awk '{print $5}' ${physioName}_extracted.txt > ${funcName}_cardio.txt
		#awk '{print $6}' ${physioName}_extracted.txt > ${funcName}_resp.txt
                nrLinesToKeep=$(echo "$(cat ${funcName}_cardio.txt|wc -l)/5"|bc)
                echo $nrLinesToKeep
		awk 'NR == 1 || NR % 5 == 0' ${funcName}_cardio.txt|head -n $nrLinesToKeep > ${funcName}_cardio_ds.txt
                nrLinesToKeep=$(echo "$(cat ${funcName}_resp.txt|wc -l)/5"|bc)
                echo $nrLinesToKeep
		awk 'NR == 1 || NR % 5 == 0' ${funcName}_resp.txt|head -n $nrLinesToKeep > ${funcName}_resp_ds.txt
		run_McRetroTS $MCRROOT Opt.Respfile=rest_resp_ds.txt Opt.Cardfile=rest_cardio_ds.txt Opt.VolTR=2.4 Opt.Nslices=43 Opt.PhysFS=100 Opt.SliceOrder=seq+z
	fi
done
