# 3D FMM<sub>ecg</sub> model overview

This repository provides a collection of functions, in the programming language R, to analyze multi-lead electrocardiogram (ECG) signals using the 3D FMM<sub>ecg</sub> model [1]. This model is built under the general assumption that the electric field of the heart is a 3-dimensional process and that the 12-lead ECG signals are the projections of that process in different directions. The 3D FMM<sub>ecg</sub> model characerizes the morphology of the five fundamental waves ($P,$ $Q,$ $R,$ $S$ and $T$) in ECG signals in terms of FMM parameters: $A$ (amplitude), $\alpha$ (phase location), $\beta$ (skewness or upward/downward peak direction) and $\omega$ (kurtosis or broadness) [2]. Indeed, there are a set of FMM parameters that are common to all the leads representing the electric field, and others that are lead-specific, representing how the signal is observed in that given direction. Moreover, the 3D FMM<sub>ecg</sub> model accurately reproduces realistic 12-lead ECG signals from healthy or pathological hearts. All of this makes the 3D FMM<sub>ecg</sub> model especially useful for the automatic diagnosis of cardiovascular diseases, patient follow-up, or decision-making on new therapies.

## How to use

Users must load functions to preprocess ECG data and fit the 3D FMM<sub>ecg</sub> model.
```
path <- getwd() # Any desired path
source(paste0(path, "/runPreprocessing_v4.1.R")) # Data preprocessing
source(paste0(path, "/FMM_ECG3D_Codes/auxMultiFMM_ECG.R")) # Data analysis
```
Preprocessing functions return preprocessed data, data segmentation in single hearthbeats and provide QRS annotations.

### Fitting example

Run the code in `fittingExample.R` which allows to fit ECG heathbeats of the patients #1 and #2 from PTB-XL database (https://physionet.org/content/ptbdb/1.0.0/).

### `fitMultiFMM_ECG` function arguments

* vDataMatrix: double matrix containing 12 columns corresponding to the standard ECG leads recorded from a heartbeat. NAs are used for an unavailable leads.
* annotation: integer indicating the QRS annotation (from 1 to the heathbeat length).
* maxIter: integer indicating the maximum number of iterations of the 3D FMM backfitting algorithm.
* parallelize: boolean. If True, a parallelized version of the fitting is performed.

## References

[1] Rueda, C., Rodríguez-Collado, A., Fernández, I., Canedo, C., Ugarte, M. D., & Larriba, Y. (2022). A unique cardiac electrocardiographic 3D model. Toward interpretable AI diagnosis. iScience, 25(12), 105617. https://doi.org/10.1016/j.isci.2022.105617

[2] Rueda, C.,  Larriba, Y., & Lamela, A. (2021). The hidden waves in the ECG uncovered revealing a sound automated interpretation method. Scientific Reports, 11, 3724. https://doi.org/10.1038/s41598-021-82520-w
