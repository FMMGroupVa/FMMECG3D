# 3D FMM<sub>ecg</sub> model codes

This repository provides a collection of functions to perform a 3D FMM<sub>ecg</sub> model [1]. The model accurately reproduces the electrocardiogram signals of any diseased or healthy heart. Furthermore, a novel algorithm accurately identifies the model parameters. This new discovery represents a revolution in electrocardiography research, solving one of the main problems in this field. It is especially useful for the automatic diagnosis of cardiovascular diseases, patient follow-up or decision-making on new therapies.

## How to use

User must load functions to preprocess ECG data and fit the 3D FMM<sub>ecg</sub> model.
```
path <- getwd() # Any desired path
source(paste0(path, "/runPreprocessing_v4.1.R")) # Data preprocessing
source(paste0(path, "/FMM_ECG3D_Codes/auxMultiFMM_ECG.R")) # Data analysis
```
Preprocessing functions returns preprocessed data, segmentation in single hearthbeats and required QRS annotations for the analysis (see **Fitting example** section).

### `fitMultiFMM_ECG` function arguments

* vDataMatrix: double matrix containing 12 columns corresponding to the standard ECG leads recorded from a heartbeat. If some leads are not availible, values shall be NAs.
* annotation: integer indicating the QRS annotation.
* maxIter: integer indicating the maximum number of iterations of the 3D backfitting algorithm.
* parallelize: boolean. When True, a parallelized version of the fitting is performed.

### Fitting example

Run the code contained in `fittingExample.R`. It allows user to fit ECG beats of the patients #1 and #2 from PTB-XL ECG database (https://physionet.org/content/ptbdb/1.0.0/).

## References

[1] Rueda, C., Rodríguez-Collado, A., Fernández, I., Canedo, C., Ugarte, M. D., & Larriba, Y. (2022). A Unique Cardiac Electrophysiological 3D Model. arXiv preprint arXiv:2202.03938.
