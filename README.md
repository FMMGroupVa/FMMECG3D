# 3D FMM<sub>ecg</sub> 

This repository provides a collection of functions, in the programming language R, to analyze multi-lead electrocardiogram (ECG) signals using the 3D FMM<sub>ecg</sub> model [1]. 

## Overview

The 3D FMM<sub>ecg</sub> model is built under the general assumption that the electric field of the heart is a 3-dimensional process and that the 12-lead ECG signals are the projections of that process in different directions. On the one hand, the 3D FMM<sub>ecg</sub> model characerizes the morphology of the five fundamental waves, $P$, $Q$, $R$, $S$ and $T$, of ECG signals in terms of FMM parameters: $A$ (amplitude), $\alpha$ (location), $\beta$ (skewness or upward/downward peak direction) and $\omega$ (kurtosis or broadness) [2]. Indeed, there are a set of FMM parameters that are common to all the leads representing the electric field, and others that are lead-specific, representing how the signal is observed in that given direction. On the other hand, the 3D FMM<sub>ecg</sub> model accurately reproduces realistic 12-lead ECG signals from healthy or pathological hearts. Moderover, the 3D FMM<sub>ecg</sub> model is especially useful for the automatic diagnosis of cardiovascular diseases, patient follow-up, or decision-making on new therapies.

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

## NORM patients analysis from PTB-XL. 3D FMM<sub>ecg</sub> Percentile Ranges Indices

PTB-XL database was analyzed using 3D FMM<sub>ecg</sub> model
|X5._P|X95._P|Mean_P|Sd_P|X5._Q|X95._Q|Mean_Q|Sd_Q|X5._R|X95._R|Mean_R|Sd_R|X5._S|X95._S|Mean_S|Sd_S|X5._T|X95._T|Mean_T|Sd_T
I_A|22,002|73,627|44,938|15,852|57,324|262,948|138,689|64,943|154,025|638,104|371,694|149,128|34,932|218,233|107,977|57,965|51,812|192,935|113,685|43,848
II_A|30,253|109,85|66,787|24,282|46,149|254,907|132,262|66,475|183,181|722,084|420,162|166,09|42,612|254,059|128,446|66,519|63,437|220,038|130,779|48,231
V1_A|13,868|59,663|33,515|14,337|61,708|300,483|159,988|76,569|173,454|757,599|435,613|181,521|67,354|372,525|196,125|93,977|19,53|134,936|63,486|38,089
V2_A|15,96|72,734|39,336|17,644|103,525|527,64|273,484|133,914|217,064|1189,753|616,442|303,316|147,713|636,997|362,287|149,064|60,272|450,34|228,606|119,734
V3_A|17,814|70,079|39,857|16,37|105,404|479,718|258,415|118,248|203,131|1061,654|573,514|270,119|166,47|644,52|373,345|149,792|84,828|464,813|247,689|117,127
V4_A|19,753|67,329|40,85|15,005|106,309|508,348|263,392|127,063|336,559|1230,165|726,805|277,21|129,452|584,903|318,72|144,369|82,262|398,297|211,929|98,907
V5_A|20,243|64|39,969|13,754|102,544|448,654|236,717|110,61|350,096|1115,739|692,33|232,764|73,632|423,253|220,787|110,034|74,012|303,489|170,084|73,399
V6_A|19,225|60,177|37,595|12,606|75,297|333,714|178,898|83,4|272,19|873,607|544,397|184,67|46,298|297,477|149,211|80,178|56,324|229,983|129,737|55,291
I_Beta|3,065|4,617|3,846|0,475|0,082|2,246|1,237|0,64|2,977|4,581|3,525|0,49|3,734|0,823|5,426|1,042|2,744|3,544|3,215|0,252
II_Beta|3,378|4,808|4,069|0,463|6,267|2,327|1,23|0,697|2,935|4,467|3,455|0,47|3,493|0,889|5,454|1,085|2,655|3,47|3,122|0,263
V1_Beta|4,854|1,37|0,115|0,812|1,858|4,816|3,31|0,91|5,413|0,344|6,078|0,389|5,976|2,511|0,666|0,748|3,583|0,255|5,06|1,11
V2_Beta|3,802|0,537|5,453|0,917|1,381|4,153|2,559|0,866|3,922|0,106|5,518|0,752|5,723|1,822|0,41|0,72|3,243|4,009|3,608|0,315
V3_Beta|3,598|6,184|4,906|0,79|0,96|3,259|1,944|0,704|3,195|6,066|4,415|0,927|5,134|1,158|0,06|0,728|3,065|3,756|3,432|0,246
V4_Beta|3,451|5,574|4,504|0,657|0,397|2,469|1,443|0,613|3,05|5,197|3,755|0,643|5,015|0,981|6,144|0,702|2,916|3,607|3,29|0,229
V5_Beta|3,316|5,161|4,235|0,567|0,075|2,178|1,209|0,619|3,026|4,67|3,56|0,503|4,709|0,884|5,882|0,757|2,789|3,524|3,195|0,236
V6_Beta|3,244|4,958|4,071|0,528|6,165|2,069|1,088|0,647|2,956|4,295|3,412|0,416|3,931|0,607|5,359|0,887|2,659|3,471|3,126|0,259
ALL_Alpha|4,417|5,152|4,847|0,224|5,506|5,637|5,609|0,041|5,745|5,876|5,765|0,039|5,862|6,005|5,889|0,048|0,954|1,766|1,336|0,254
ALL_Omega|0,067|0,219|0,135|0,047|0,024|0,064|0,039|0,016|0,026|0,044|0,035|0,008|0,019|0,051|0,031|0,01|0,124|0,288|0,181|0,052


## References

[1] Rueda, C., Rodríguez-Collado, A., Fernández, I., Canedo, C., Ugarte, M. D., & Larriba, Y. (2022). A unique cardiac electrocardiographic 3D model. Toward interpretable AI diagnosis. iScience, 25(12), 105617. https://doi.org/10.1016/j.isci.2022.105617

[2] Rueda, C.,  Larriba, Y., & Lamela, A. (2021). The hidden waves in the ECG uncovered revealing a sound automated interpretation method. Scientific Reports, 11, 3724. https://doi.org/10.1038/s41598-021-82520-w
