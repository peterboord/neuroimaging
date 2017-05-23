#SFLRCM
for x in $(ls mni*[SFLRCMl]);do 
outFile="$(echo $x|sed 's/mniCoordMm/mniCoordVox/')"
fcmd="fslmaths $MNI1mm -mul 0 -add 1 -roi $(echo $(echo $(cat $x) | std2imgcoord -img $MNI1mm -std $MNI1mm -vox -)|sed 's/ / 1 /g') 1 0 1 -bin $outFile"
echo $cmd
eval $cmd
fslmaths $outFile -s 1 $outFile
done

