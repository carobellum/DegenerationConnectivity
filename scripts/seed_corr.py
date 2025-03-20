#!/Users/callithrix/Documents/Projects/Degeneration/code/cerebellar_degeneration/env python
# -*- coding: utf-8 -*-
"""
Calculate correlation between seed resting-state timecourses
"""


from pathlib import Path
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import util as util
import nitools as nt
import itertools
from scipy import stats
from copy import deepcopy


base_dir = '/Users/callithrix/Documents/Projects/Degeneration/derivatives/'

def prepare_atlas_ts(regions=['M3', 'D2']):
    """
    Prepare atlas time series for selected regions. Imports the ts_atlas_decon.txt file for each subject and session,
    extracts the timecourses for the selected regions, and saves the timecourses in separate txt files.

    Args:
        regions (list): List of region labels to extract timecourses for.


    """
    regions = [f'{region}L' for region in regions] + [f'{region}R' for region in regions]
    
    ts_dir = sorted(glob.glob(f'{base_dir}/sub-*/ses-*/func/ts'))
    _, _, labels = nt.read_lut(f'/Users/callithrix/code/Python/cerebellar_atlases/Nettekoven_2023/atl-NettekovenSym32.lut')
    region_indices = [labels[1:].index(region) for region in regions]
    for scan_no, dir in enumerate(ts_dir):
        # if file f'{dir}/ts_atlas_decon.txt' exists
        if Path(f'{dir}/cereb_atlas/ts_atlas_decon.txt').exists():
            t = np.loadtxt(f'{dir}/cereb_atlas/ts_atlas_decon.txt')
            for r,region_indx in enumerate(region_indices):
                np.savetxt(f'{dir}/ts_cereb-{regions[r][:-1]}-{"left" if "L" in regions[r] else "right"}.txt', t[:, region_indx])


def import_ts():
    """
    Imports the timeseries of each ROI for each subject and session and concatenates them into a 3D array.
    Needs to be run after preprare_atlas_ts(), since preprare_atlas_ts() extracts the cerebellar atlas timecourses 
    and saves them separately in the ts directory of each subject and session.
    """
    ts_dir = sorted(glob.glob(f'{base_dir}/sub-*/ses-*/func/ts'))

    for scan_no, dir in enumerate(ts_dir):
        ts_files = glob.glob(f'{dir}/*.txt')
        for roi_no, file in enumerate(ts_files):
            t = np.loadtxt(file)
            if scan_no == 0 and roi_no == 0:
                # timepoints, rois, n
                ts = np.zeros((t.shape[0], len(ts_files), len(ts_dir)))
            ts[:, roi_no, scan_no] = t
    rois = ['-'.join(file.split('/ts_')[-1].split('.txt')[0].split('_'))
            for file in ts_files]
    subjects = [
        'sub-{}'.format(dir.split('/sub-')[1].split('/')[0]) for dir in ts_dir]
    sessions = [
        'ses-{}'.format(dir.split('/ses-')[1].split('/')[0]) for dir in ts_dir]
    return ts, rois, subjects, sessions


def correlate_ts(ts):
    """
    Calculate the correlation coefficients between time series data.

    Parameters:
    ts (ndarray): The input time series data. It should have shape (num_samples, num_features, num_timepoints).

    Returns:
    ndarray: The correlation coefficients between the time series data. It has shape (num_features, num_features, num_timepoints).
    """

    for n in np.arange(ts.shape[2]):
        timecourses = ts[:, :, n]
        # Standardize the time series for easier calculation
        timecourses = util.zstandarize_ts(timecourses)
        coef = timecourses.T @ timecourses / timecourses.shape[0]
        # divide
        if n == 0:
            coefs = np.zeros((ts.shape[1], ts.shape[1], ts.shape[2]))
        coefs[:, :, n] = coef

    return coefs


def get_conn(coefs, rois):
    """
    Calculate the connectivity matrix from the correlation coefficients.

    Parameters:
    coefs (ndarray): The correlation coefficients matrix of shape (n, n, m),
                     where n is the number of regions of interest (ROIs) and m is the number of samples.
    rois (list): The list of ROIs.

    Returns:
    ndarray: The connectivity matrix of shape (m, k), where k is the number of unique pairs of ROIs.

    """
    # get off-diagonal elements
    pairs = ['_'.join(pair) for pair in itertools.combinations(
        rois, 2) if not pair[0] == pair[1]]
    for i in np.arange(coefs.shape[2]):
        np.triu(coefs[:, :, i])
        coef = coefs[:, :, i]
        coef = coef[np.triu_indices(coef.shape[0], k=1)]
        if i == 0:
            conn = np.zeros((coefs.shape[2], len(pairs)))
        conn[i, :] = coef

    return conn


def get_data(subjects, sessions, rois, conn):
    """
    Retrieves and processes data for connectivity analysis.

    Args:
        subjects (list): List of subject IDs.
        sessions (list): List of session IDs.
        rois (list): List of regions of interest.
        conn (numpy.ndarray): Connectivity data.

    Returns:
        tuple: A tuple containing two pandas DataFrames:
            - data_long: Long-format DataFrame with subject-level connectivity data.
            - data_wide: Wide-format DataFrame with subject-level connectivity data, demographic data, behavioural data and additional information.

    """

    # get connectivity pair names
    pairs = ['_'.join(pair) for pair in itertools.combinations(
        rois, 2) if not pair[0] == pair[1]]

    # make dataframe
    d = {'subject': subjects, 'session': sessions}
    connectivity = pd.concat(
        [pd.DataFrame(d), pd.DataFrame(conn, columns=pairs)], axis=1)
    connectivity.session = connectivity.session.astype(
        "category").cat.reorder_categories(['ses-pre', 'ses-post'])

    # import patient data
    info = pd.read_csv(f'{base_dir}/participants.tsv', sep='\t')

    # import training data
    training = pd.read_csv(f'{base_dir}/phenotype/trainingdata.tsv', sep='\t')
    data = pd.merge(connectivity, info, left_on='subject',
                    right_on='participant_id', how='inner').drop(columns=['participant_id'])
    cols_to_use = training.columns.difference(data.columns)
    data_wide = pd.merge(data, training[cols_to_use], left_on='subject',
                         right_on='participant_id', how='inner')
    # drop superfluous participant index and ses-pre for sub-57, which has no measurements (all nan)
    data_wide = data_wide.drop(columns=['participant_id']).dropna(
        subset=connectivity.columns[2:], how='all')

    # add proprioceptive and feedback condition indicator
    cond_inf = data["condition"].str.split("F", expand=True)
    cond_inf.columns = ["vision", "feedback"]
    cond_inf.vision[cond_inf.vision == "VM"] = True
    cond_inf.vision[cond_inf.vision == "Prop"] = False
    cond_inf.feedback[cond_inf.feedback == "B"] = True
    cond_inf.feedback[cond_inf.feedback != True] = False

    data_wide = data_wide.join(cond_inf)
    data_long = pd.melt(data_wide, id_vars=['subject', 'session', 'age', 'sex', 'group', 'condition',
                                            'BA6_L_Post_Int', 'BA6_L_Pre_Int', 'Diff_CuneusR', 'Diff_PMd_L_Interac',
                                            'Diff_PMd_R_Pat_patMask', 'Diff_RJPE_A10', 'Diff_RJPE_A50',
                                            'PAVG_A10_Post', 'PAVG_A10_Pre', 'PAVG_A25_Post', 'PAVG_A25_Pre',
                                            'PAVG_A50_Post', 'PAVG_A50_Pre', 'PSTD_A10_Post', 'PSTD_A10_Pre',
                                            'PSTD_A25_Post', 'PSTD_A25_Pre', 'PSTD_A50_Post', 'PSTD_A50_Pre',
                                            'Post_BA6R_patMask', 'Post_CuneusR', 'Pre_BA6R_patMask', 'Pre_CuneusR',
                                            'SARA', 'TIV', 'vision', 'feedback'],
                        value_vars=pairs, value_name='connectivity', var_name='regions')

    return data_long, data_wide

def regions_main_analysis():
    regions_contra = [
                    'ppc-right_cereb-M3-left',
                    'ppc-left_cereb-M3-right',
                    'pmd-left_cereb-M3-right',
                    'pmd-right_cereb-M3-left',
                    'm1-hand-left_cereb-M3-right',
                    'm1-hand-right_cereb-M3-left',
                    'dlpfc-right_cereb-D2-left',
                    
                    ]

    regions_ipsi = [
                    'ppc-right_cereb-M3-right',
                    'ppc-left_cereb-M3-left',
                    'pmd-right_cereb-M3-right',
                    'pmd-left_cereb-M3-left',
                    'm1-hand-right_cereb-M3-right',
                    'm1-hand-left_cereb-M3-left',
                    'dlpfc-right_cereb-D2-right',
                    
                    ]

    regions_cortico = [
                    'ppc-left_ppc-right',
                    'pmd-left_pmd-right',
                    'm1-hand-left_m1-hand-right',
                    ]

    regions_cereb = [
                    'cereb-D2-left_cereb-D2-right',
                    'cereb-M3-left_cereb-M3-right',
                    ]
    regions_control = [
                    'dlfpc-right_cereb-S1-left',
                    'dlfpc-right_cereb-S1-right',                    
                    ]
    
    regions_selected = [*regions_contra, *
                        regions_ipsi, *regions_cortico, *regions_cereb, *regions_control]
    return regions_selected, [regions_contra, regions_ipsi, regions_cortico, regions_cereb, regions_control]

def regions_control_analysis():
    regions_control = [
                    'dlfpc-right_cereb-S1-left',
                    'dlfpc-right_cereb-S1-right',                    
                    ]
    
    regions_selected = [*regions_control]
    return regions_selected, [regions_control]

def select_regions(data_long, data_wide, control_regions=False):
    """
    Selects specific regions from the input data based on predefined criteria.

    Args:
        data_long (DataFrame): The long-format data containing region information.
        data_wide (DataFrame): The wide-format data containing region information.

    Returns:
        Tuple: A tuple containing the following:
            - data_long (DataFrame): The updated long-format data with selected regions.
            - data_wide (DataFrame): The updated wide-format data with selected regions.
            - regions_contra (list): A list of selected contra-lateral regions.
            - regions_ipsi (list): A list of selected ipsi-lateral regions.
            - regions_cortico (list): A list of selected cortico-cortical regions.
            - regions_cereb (list): A list of selected cerebello-cerebellar regions.
    """

    # Rename region names to naming convention: < cortical > - < hemishpere > _ < cerebellar > - < hemisphere >
    regions = data_long['regions'].unique()

    # Reorder regions so that cortical region comes first and left comes before right
    cortical = ['m1', 'pmd', 'ppc', 'dlpfc',
                'cereb']
    hemispheres = ['left', 'right']

    # Reordering
    for r, region in enumerate(regions):
        region_elements = region.split('_')

        # reorder regions
        position1 = cortical.index(region_elements[0].split('-')[0])
        position2 = cortical.index(region_elements[1].split('-')[0])

        if position1 > position2:
            new_region = '_'.join([region_elements[1], region_elements[0]])
        elif position1 == position2 and hemispheres.index(region_elements[0].split(
                '-')[-1]) > hemispheres.index(region_elements[1].split('-')[-1]):
            new_region = '_'.join([region_elements[1], region_elements[0]])
        else:
            new_region = '_'.join([region_elements[0], region_elements[1]])
        regions[r] = new_region

    replace_regions = dict(zip(data_long['regions'].unique(), regions))

    data_wide = data_wide.rename(replace_regions, axis='columns')
    data_long['regions'] = data_long['regions'].replace(replace_regions)

    if control_regions==True:
        regions_selected, region_list = regions_control_analysis()
    else:
        regions_selected, region_list = regions_main_analysis()
    

    data_long = data_long[data_long['regions'].isin(regions_selected)]

    selected_columns = [*data_long.columns.drop('regions').drop('connectivity'),
                        *regions_selected]
    data_wide = data_wide[selected_columns]

    return data_long, data_wide, region_list


def select_atlas_regions(ts, rois, selected_rois):
    pass
    selected_rois = [roi for roi in selected_rois if roi in rois]
    selected_rois_idx = [rois.index(roi) for roi in selected_rois]
    ts = ts[:, selected_rois_idx, :]
    rois = [rois[i] for i in selected_rois_idx]
    return ts, rois


if __name__ == '__main__':
    prepare_atlas_ts(regions=['S1'])
    ts, rois, subjects, sessions = import_ts()
    coefs = correlate_ts(ts)
    conn = get_conn(coefs, rois)
    data_long, data_wide = get_data(subjects, sessions, rois, conn)
    data_long, data_wide, region_list = select_regions(
        data_long, data_wide)
    regions_contra, regions_ipsi, regions_cortico, regions_cereb, regions_control = region_list
    data_long.to_csv(f'data/connectivity.tsv', sep='\t')
    data_wide.to_csv(f'data/connectivity_wide.tsv', sep='\t')


    data_long, data_wide, region_list = select_regions(
    data_long, data_wide, control_regions=True)
    # regions_contra, regions_ipsi, regions_cortico, regions_cereb, regions_control = region_list
    data_long.to_csv(f'data/connectivity_control.tsv', sep='\t')
    data_wide.to_csv(f'data/connectivity_wide_control.tsv', sep='\t')

    