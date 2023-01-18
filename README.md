# 3D FMM<sub>ecg</sub> 

This repository provides a collection of functions to analyse electrocardiogram (ECG) signals using the 3D FMM<sub>ecg</sub> model [1]. It is based on the novel  Frequency Modulated Möbius (FMM) approach recently developed by our group (http://www.eio.uva.es/the-fmm-project/). The code is developed in the programming language R and requires of the R package FMM [2,3].

## Overview

3D FMM<sub>ecg</sub> is an interpretable, flexible, universal, and freely available approach for the automatic interpretation of ECG. 
The 3D FMM<sub>ecg</sub> model is built under the general assumption that the electric field of the heart is a 3-dimensional process and that the 12-lead ECG signals are the projections of that process in different directions. On the one hand, the 3D FMM<sub>ecg</sub> model characerizes the morphology of the five fundamental waves, $P$, $Q$, $R$, $S$ and $T$, of ECG signals in terms of FMM parameters: $A$ (amplitude), $\alpha$ (location), $\beta$ (skewness or upward/downward peak direction) and $\omega$ (kurtosis or width), see Figure 1. Indeed, there are a set of FMM parameters ( $\alpha$, $\omega$ ) that are common to all the leads representing the electric field, and others ( $A$, $\beta$ ) that are lead-specific, representing how the signal is observed in that given direction. On the other hand, the 3D FMM<sub>ecg</sub> model accurately reproduces realistic multi-lead ECG signals from healthy or pathological hearts. This property allows to elucidate how any algorithm operates, detecting possible shortcoming problems and identifying the causes of its decisions. Moreover, 3D FMM<sub>ecg</sub> model works regardless the recording device,
the number of leads, the length of data, or the differences between datasets label distributions. 

<p align="center">
  <img src="https://user-images.githubusercontent.com/117477025/213187046-2fb6652a-53ed-4f6d-8e1e-ae7f91571e66.jpg" width="450" height="300" alt>
</p>

Figure 1: Typical ECG heartbeat. Wave decomposition and morphological FMM parameter description. $\alpha$, $\beta \in [0,2 \pi)$, $\omega \in (0,1]$  and $A \in \mathcal{R}^+.$

## How to use
3D FMM<sub>ecg</sub> is designed to analyse multi-lead ECG fragments of any length. Data are preprocessed to  achieve reliable ECG fragments and to divide them beat by beat. Then, for each heartbeat, 3D FMM<sub>ecg</sub> provides FMM parameter estimates. For each wave, the series of parameter values, corresponding to consecutive beats, which can be summarized to get average patterns as well as the changes in the patterns over time.

### Data Preprocessing

3D FMM<sub>ecg</sub> code incorporates a standard preprocessing for single or multi-lead ECG data including baseline correction, QRS detection and ECG segmentation.

Users must load functions to preprocess ECG data and run `givePreprocessing_app` function.

#### `givePreprocessing_app` function arguments

* dataIn: double matrix of 12 columns with raw ECG fragment data across leads. Lead names must be provided as header. NAs are used for an unavailable leads.
* freqHz: integer indicating indicating the sampling rate in Hertz (Hz).

This function returns the preprocessed data, QRS annotations, ECG segmentation and lead-specific unreliable heartbeats to be discarded from the analysis.

### Data Analysis

Prepocessed ECG fragment is divided into beats to be subsequently analysed using 3D FMM<sub>ecg</sub> model. For computational efficiently, only the leads: I, II, V1, V2, V3, V4, V5 and V6 are considered for the analysis. The rest are linear combinations of I and II.

Users must load functions to analyse ECG data and run `fitMultiFMM_ECG` function.

#### `fitMultiFMM_ECG` function arguments

* vDataMatrix: double matrix of 12 columns with heartbeat preprocessed ECG data across leads. NAs are used for an unavailable leads.
* annotation: integer indicating the QRS annotation (from 1 to the heathbeat length).

This function returns FMM wave parameter estimates across the eight leads considered.

### Fitting example

Run the code in `fittingExample.R` to analyse using 3D FMM<sub>ecg</sub> patient #1 from PTB-XL database (https://physionet.org/content/ptb-xl/1.0.3/) [4].

## NORM patient analysis from PTB-XL. Percentile Ranges for 3D FMM<sub>ecg</sub> Indices

PTB-XL is a large dataset of 21837 clinical 12-lead ECGs of 10 second length annotated by two cardiologists with diagnostic labels, based on SCP-ECG statements, and the likelihood information for the  statements [4]. PTB-XL database has been analysed using 3D FMM<sub>ecg</sub> in [1] for patients with likelihood $\geq 80$. In particular, we analysed 9055 patients from PTB-XL labelled as NORM, i.e. with normal ECGs. 

The normal percentile ranges ($5th}$, $95th$) of several related 3D FMM<sub>ecg</sub> indices calculated from NORM patients are specially useful for identifyng noisy and/or pathological ECG patters, as those with values out such ranges. The FMM-based incdices for which these ranges were computed are: 

* FMM parameterS: $A, \beta$ are lead-specific. $\alpha, \omega$ are equal across leads.
* $Var_J$: measure of the relative relevance of wave $J$, for $J= P, Q, R, S, T$. For a given lead, it is calculated as the variability the wave $J$ explains, see [5] for details.
* $R^2:$ measure of the accuracy of the model across leads, see [5] for details.
* RR: duration of R-R interval from QRS annotations in milisecons (ms).
* disPQ, disQS, and disQT: difference in ms between $\alpha_P$ and $\alpha_Q;$, $\alpha_Q$ and $\alpha_S$; and $\alpha_S$ and $\alpha_T$, respectively.

Normal percentile ranges for the median (Me) and coeficcient of variation (Cv) of these indices across NORM patients in PTB-XL are given in Tables 1-5 and 6-10, respectively. Cv for angular parameters was defined in the Suplementary Material of [5]. 


||$5th_P$|$95th_P$|$Mean_P$|$Sd_P$|$5th_P$|$95th_P$|$Mean_P$|$Sd_P$|$5th_P$|$95th_P$|$Mean_P$|$Sd_P$|$5th_P$|$95th_P$|$Mean_P$|$Sd_P$|$5th_P$|$95th_P$|$Mean_P$|$Sd_P$|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---
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

Table 1: $5th$, $95th$, Mean and standard deviation (Sd) for the Me of FMM parameters. Subindes in header denote waves. ALL is used to denote a common estimation across leads. 


## References

[1] Rueda, C., Rodríguez-Collado, A., Fernández, I., Canedo, C., Ugarte, M. D., & Larriba, Y. (2022). A unique cardiac electrocardiographic 3D model. Toward interpretable AI diagnosis. iScience, 25(12), 105617. https://doi.org/10.1016/j.isci.2022.105617

[2] Fernández, I., Rodríguez-Collado, A., Larriba, Y., Lamela, A., Canedo, C., and Rueda, C. (2021). FMM: rhythmic patterns modeling by FMM models. R package version 0.3.0.

[3] Fernández, I., Rodríguez-Collado, A., Larriba, Y., Lamela, A., Canedo, C., and Rueda, C. (2022). FMM: an R package for modeling rhythmic patterns in oscillatory systems. R Journal. 14, 361–380.

[4] Wagner, P., Strodthoff, N., Bousseljot, RD. et al. (2020) PTB-XL, a large publicly available electrocardiography dataset. Scientific Data, 7, 154. https://doi.org/10.1038/s41597-020-0495-6

[5] Rueda, C.,  Larriba, Y., & Lamela, A. (2021). The hidden waves in the ECG uncovered revealing a sound automated interpretation method. Scientific Reports, 11, 3724. https://doi.org/10.1038/s41598-021-82520-w

