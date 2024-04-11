#!/bin/bash
# ------------------------------------------------------------------------------
# Script name:  seed_ts.sh
#
# Description:
#               Script to extract seed timecourses from cleaned functional data
#
# Author:       Caroline Nettekoven, 2022
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -d ~/Documents/Projects/Cerebellar_Degeneration/code/cerebellar_degeneration ]; then
SCRIPTDIR=~/Documents/Projects/Cerebellar_Degeneration/code/cerebellar_degeneration
WORKDIR=~/Documents/Projects/Cerebellar_Degeneration/derivatives/
hard_drive_workdir=/Volumes/nettekoven_data/Projects/Cerebellar_Degeneration/derivatives
# rois=('m1_hand_right' 'm1_hand_left' 'cereb_hand_left_hardw' 'cereb_hand_right_hardw')
rois=('pmd_left' 'pmd_right' 'sma_left' 'sma_right' 'cereb_hand_left_mdtb' 'cereb_hand_right_mdtb')
rois=('ppc_left' 'ppc_right' 'dlpfc_right' 'cereb_hand_left_mdtb' 'cereb_hand_right_mdtb' 'crus1_right' 'crus2_right' 'lobule6_right' 'crus1_left' 'crus2_left' 'lobule6_left')
sessions=('pre' 'post')
elif [ -d /home/fs0/cnette/scratch/Cerebellar_Degeneration/code/cerebellar_degeneration ]; then
SCRIPTDIR=/home/fs0/cnette/scratch/Cerebellar_Degeneration/code/cerebellar_degeneration
WORKDIR=/home/fs0/cnette/scratch/Cerebellar_Degeneration/derivatives/
datadir=/vols/Data/ping/caro/backups/degeneration/derivatives
template_dir=/home/fs0/cnette/scratch/Cerebellar_Degeneration/derivatives//template
roi_dir=/vols/Data/ping/caro/backups/degeneration/derivatives/rois
# rois='m1_hand_right m1_hand_left'
# rois=' cereb_hand_left_mdtb  cereb_hand_right_mdtb'
# rois=' pmd_left  pmd_right'
# rois=' ppc_left  ppc_right dlpfc_right'
rois=' atlas_decon'
subjects=`ls -d ${datadir}/sub-*`
sessions='pre post'
else
echo "Cerebellar_Degeneration working directory not found."
fi
# ----
cd ${SCRIPTDIR}

# ------ Extract timeseries from seeds ------
for i in $subjects; do
    subject=${i#*derivatives/}
    for session in $sessions; do
        
    
        for roi in $rois; do

        if [ ! -f $i/ses-${session}/func/ts/ts_${roi}.txt ] && [ -f $i/ses-${session}/func/rois/${roi}.nii.gz ] ; then
        echo Extracting timecourse of ${subject} ${session} ${roi}...

        # echo \
        # fsl_sub \
        # fslmeants \
        # -i $i/ses-${session}/func/filtered_func_data_clean \
        # -m $i/ses-${session}/func/rois/${roi} \
        # -o $i/ses-${session}/func/ts/ts_${roi}.txt

        # echo \
        fsl_sub \
        fslmeants \
        -i $i/ses-${session}/func/filtered_func_data_clean \
        --label=$i/ses-${session}/func/rois/atlas_decon \
        -o $i/ses-${session}/func/ts/ts_${roi}.txt


        fi

        done
    done
done

