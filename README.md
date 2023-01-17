# 3D FMM<sub>ecg</sub> model 

This repository provides a collection of functions to perform a 3D FMM<sub>ecg</sub> model [1]. The mathematical formulation of this model accurately reproduces and parametrizes ECG signals of any diseased or healthy heart. The 3D FMM<sub>ecg</sub> model is especially useful for the automatic diagnosis of cardiovascular diseases, patient follow-up or decision-making on new therapies.

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
