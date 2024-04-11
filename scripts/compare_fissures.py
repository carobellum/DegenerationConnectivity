#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare distances between cerebellar fissures
"""


from pathlib import Path
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import Functional_Fusion.util as util
import itertools
from scipy import stats
import nibabel as nb
import math
from scipy import stats


base_dir = '/Users/callithrix/Documents/Projects/Cerebellar_Degeneration/derivatives/'


def import_fissures(fissure=1, space='degeneration_40'):
    """This function is used to import fissure data from a set of NIFTI files and return the fissures data and corresponding subject names.
    Space is either degeneration_40 or mni_ants
    """
    files = sorted(
        glob.glob(f'{base_dir}/sub-*/ses-pre/anat/rois/fissure-{fissure}_{space}.nii.gz'))

    for f, file in enumerate(files):
        img = nb.load(file)
        if f == 0:
            fissures = []
        fissures.append(img)

    subjects = [
        'sub-{}'.format(file.split('/sub-')[1].split('/')[0]) for file in files]

    return fissures, subjects


def import_mask(mask='vermal', space='degeneration_40'):
    mask_file = base_dir + f'/rois/{mask}_{space}.nii.gz'
    mask = nb.load(mask_file)
    return mask


def get_coord(fissure, mask=None):
    idx = fissure.get_fdata() > 0  # find voxels that are not zeroes
    if mask is not None:
        idx = np.logical_and(idx, mask.get_fdata())
    idx = np.where(idx)
    ijk = np.vstack(idx).T  # list of arrays to (voxels, 3) array
    xyz = nb.affines.apply_affine(fissure.affine, ijk)  # get mm coords
    return xyz


def calculate_distance(point_A, point_B):
    """
    Calculates the Euclidean distance between two points in a three-dimensional space.

    """
    distance = math.sqrt((point_A[0] - point_B[0])**2 +
                         (point_A[1] - point_B[1])**2 + (point_A[2] - point_B[2])**2)

    return distance


def fissure_distance(fissure_A, fissure_B, mask):
    """
    Calculates the average minimum Euclidean distance between two fissures.

    Parameters:
    fissure_A (nibabel image object): A nibabel image object representing a fissure.
    fissure_B (nibabel image object): A nibabel image object representing another fissure.
    mask (nibabel image object): A binary mask used for calculation of fissure distance.

    Returns:
    float: The average minimum Euclidean distance between fissure_A and fissure_B.
    """
    coord_A = get_coord(fissure_A, mask=mask)
    coord_B = get_coord(fissure_B, mask=mask)

    # Get the number of points in each fissure
    points_A = coord_A.shape[0]
    points_B = coord_B.shape[0]

    if points_A == 0 or points_B == 0:
        average_distance = np.nan
    else:
        # Loop through the points in fissure A
        D = []
        for i in range(points_A):
            point = coord_A[i, :]
            for j in range(points_B):
                if j == 0:
                    distance = np.zeros(points_B)
                distance[j] = calculate_distance(point, coord_B[j, :])
            D.append(distance.min())
        average_distance = sum(D) / len(D)

    return average_distance


def get_distances(subjects, fissures, mask):
    """
    Calculates the average distances between fissures for each subject.

    Parameters:
    subjects (list): A list of subjects.
    fissures (list): A list of fissures for each subject.
    mask (ndarray): A binary mask used for calculation of fissure distance.

    Returns:
    list: A list of average distances between fissures for each subject.

    Example:
    >>> subjects = ['subject1', 'subject2', 'subject3']
    >>> fissures = [fissure1, fissure2, fissure3]
    >>> mask = np.array([[0, 0, 1], [0, 1, 1], [1, 1, 1]])
    >>> get_distances(subjects, fissures, mask)
    [0.5, 1.0, 0.75]
    """

    Distance = []
    for i in range(len(fissures)):
        print(f'Subject {subjects[i]}')
        distances = []
        for j in range(len(fissures)):
            # ignore distances between the same subject
            if i == j:
                distance = np.nan
            else:
                distance = fissure_distance(fissures[i], fissures[j], mask)
            distances.append(distance)
        distances = [dist for dist in distances if ~np.isnan(dist)]
        if np.all(np.isnan(distances)):
            average_distance = np.nan
        else:
            average_distance = sum(distances) / len(distances)
        Distance.append(average_distance)

    return Distance


def make_dataframe(fissure=1):
    if fissure == 1:
        mask_names = ['vermal', 'cerebellar_mask']
    elif fissure == 2:
        mask_names = ['cerebellar_mask']

    data = []
    for mask_name in mask_names:
        # for space in ['pc', 'p', 'c', 'pc_cerebellum', 'p_cerebellum', 'c_cerebellum']:
        for space in ['pc', 'p', 'c']:
            fissures, subjects = import_fissures(
                fissure=fissure, space=space)
            mask = import_mask(mask=mask_name, space=space)
            distances = get_distances(subjects, fissures, mask)

            # make dataframe
            d = {'subject': subjects,
                 'fissure': [fissure] * len(subjects),
                 'mask': [mask_name] * len(subjects),
                 'space': [space] * len(subjects),
                 'distances': distances}
            data.append(pd.DataFrame(d))
    Data = pd.concat(data)

    return Data


if __name__ == '__main__':

    D1 = make_dataframe(fissure=1)
    D2 = make_dataframe(fissure=2)
    D = pd.concat([D1, D2])
    D.to_csv(f'{base_dir}/data/fissures_new.csv', index=False, sep='\t')
    pass

