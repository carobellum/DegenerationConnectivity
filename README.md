# DegenerationConnectivity
Code for analysing functional connectivity of cerebellar degeneration patients for manuscript Nettekoven et al., 2024.

## Dependencies:

see ```requirements.txt```

## Dataset information
The data files available in this repository were derived from MRI scans of 40 patients diagnosed with pure cerebellar cortical degeneration and 40 age and sex-matched neurologically healthy individuals. All individuals participated a five-day motor training and underwent a functional and structural MRI scan on the days before and after training. 

Functional data is missing for (subject, timepoint):
- sub-57, post
- sub-70, pre
- sub-70, post

Structural data was acquired from all subjects.

## Template information
A study-specific template is included into this repository under template/DeCon_template.nii.gz

| File                                  | Description                                                |
| ------------------------------------- | ---------------------------------------------------------- |
| DeCon_template.nii.gz                 | DeCon template image                                       |
| DeCon_template_mask.nii.gz            | DeCon template whole-brain mask                            |
| DeCon_template_cerebellar_mask.nii.gz | DeCon template cerebellar mask                             |
| DeCon_to_MNI_warp.nii.gz              | DeCon template non-linear warp to MNI152NLin6Asym template |
| MNI_to_DeCon_warp.nii.gz              | MNI152NLin6Asym template non-linear warp to DeCon template |


## Notebooks / Code to replicate different sections of the paper

### Study demographics
Study demographics, demographics for template generation and validation sample (Supplementary Table 1), and demographics for FIX training datasets
```demographics.ipynb```

### Study-specific template
Fissure distances in template, SUIT and MNI space was calculated using:
```compare_fissures.py```

Fissure overlap was plotted (Fig 1B-D) and compared using:
```stats_template.ipynb```

FIX performance with template registrations and standard registrations was plotted (Fig 1E) and evaluated using:
```stats_fix.ipynb```

### Connectivity

#### Extracting ROI timecourses

ROI timecourses were extracted from functional data in native space using:
```seed_ts.sh```

#### Correlating ROI timecourses

Functional connectivity between timecourses was calculated using:
```seed_corr.py```

#### Statistical analysis of functional connectivity
Data was loaded from dataframes, normalized and brought into the correct shape for analysis using:
```data_connectivity.R```

Baseline connectivity differences were plotted (Fig 3A & 3B) using:
```plots_connectivity.ipynb```

Connectivity change was plotted (Fig 4, 5, 6 & 7) using:
```plots_connectivity_change.ipynb```

Statistical tests on functional connectivity were calculated using:
```stats.R```

Model assumptions were tested using:
```model_assumptions.R```

Connectivity results were plotted using:
```plots_connectivity.R``` 

