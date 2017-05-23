dirNames=$(wildcard RC4[1-2][0-9][0-9])
targets=$(dirNames:%=/projects2/udall/task/%/session1/T1_to_MNI_deformed.nii.gz)




# Set the default number of threads used by ANTS to 1 to parallelize on 
# gridengine
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1

#unset ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
#targets=/projects2/udall/task/RC4106/session1/T1_to_MNI_deformed.nii.gz




all: $(targets)

%/T1_to_MNI_deformed.nii.gz:
	${ANTSPATH}antsIntroduction.sh -d 3 -i $*/T1_brain.nii.gz -m 30x90x20 -o $*/T1_to_ANTStemplate_ -s CC -r /projects2/udall/ANTS/udall_template.nii.gz -t GR

write/seeds/zstat%_write_Alphabetical-order.feat: write/seeds/write_ssmooth_Alphabetical-order.nii.gz write/seeds/thresh_zstat%_write_Alphabetical-order.txt $(FEAT_TASK_TEMPLATE) /usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz $(STANDARD_DIR)/NIHtoMNIWarp.nii.gz $(STANDARD_DIR)/NIHtoMNIAffine.txt xfm_dir/T1_to_nih_deformed.nii.gz xfm_dir/write_to_nih_epireg_ants.nii.gz
        rm -rf $@ ;\
        nvols=`fslval $(word 1,$^) dim4` ;\
        sed -e 's/SUBJECT/$(SUBJECT)/g' -e 's/SESSION/$(SESSION)/g' -e 's/TASK/'write'/g' -e 's/RUN/'Alphabetical-order'/g' -e "s/NVOLS/$$nvols/g" -e "s/COMPONENT/`basename $(word 2,$^)`/g" -e "s/OUTDIR/`basename $@`/g" $(FEAT_TASK_TEMPLATE) > `dirname $@`/`basename $@`.fsf ;\
        export FSLPARALLEL=false; feat `dirname $@`/`basename $@`.fsf ;\
        export ANTSPATH=$(ANTSpath) ;\
        $(ANTSpath)/WarpImageMultiTransform 3 $@/stats/cope1.nii.gz $@/stats/cope1_in_std.nii.gz -R $(word 4,$^) $(word 5,$^) $(word 6,$^) xfm_dir/`basename $(word 7,$^) _deformed.nii.gz`_Warp.nii.gz xfm_dir/`basename $(word 7,$^) _deformed.nii.gz`_Affine.txt xfm_dir/write_to_T1_ras.txt
