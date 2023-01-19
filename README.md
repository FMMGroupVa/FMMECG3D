# 3D FMM<sub>ecg</sub> 

This repository provides a collection of functions to analyse one or multi-lead electrocardiogram (ECG) signals using the 3D FMM<sub>ecg</sub> model [1]. It is based on an efficiently mathematical formulation, based on Frequency Modulated Möbius (FMM), recently developed by our group (http://www.eio.uva.es/the-fmm-project/). 
The code is developed in the programming language R and requires of the R package FMM [2,3].

## Overview

3D FMM<sub>ecg</sub> is an interpretable, flexible, universal, and freely available approach for the automatic interpretation of ECG. 
The 3D FMM<sub>ecg</sub> model is built under the general assumption that the electric field of the heart is a 3-dimensional process and that the 12-lead ECG signals are the projections of that process in different directions. On the one hand, the 3D FMM<sub>ecg</sub> model characerizes the morphology of the five fundamental  waves, labelled as  $P$, $Q$, $R$, $S$ and $T$  and located in order,  of ECG signals in terms of the FMM parameters: $A$ (amplitude), $\alpha$ (location), $\beta$ (skewness or upward/downward peak direction) and $\omega$ (kurtosis or width), see Figure 1. Indeed, there are a set of FMM parameters ( $\alpha$, $\omega$ ) that are common to all the leads representing the electric field, and others ( $A$, $\beta$ ) that are lead-specific, representing how the signal is observed in that given direction. On the other hand, the 3D FMM<sub>ecg</sub> model accurately reproduces realistic multi-lead ECG signals from healthy or pathological hearts. This property allows to elucidate how any algorithm operates, detecting possible shortcoming problems and identifying the causes of its decisions. Moreover, 3D FMM<sub>ecg</sub> model works regardless the recording device, the number of leads, the length of data, or the differences between datasets label distributions. 

<p align="center">
  <img src="https://user-images.githubusercontent.com/117477025/213187046-2fb6652a-53ed-4f6d-8e1e-ae7f91571e66.jpg" width="450" height="300" alt>
</p>

Figure 1: Typical ECG heartbeat. Wave decomposition and morphological wave-specific FMM parameter description. $\alpha$, $\beta \in [0,2 \pi)$, $\omega \in (0,1]$  and $A \in \mathcal{R}^+.$

## How to use
3D FMM<sub>ecg</sub> is designed to analyse one or multi-lead ECG fragments of any length. Data are preprocessed to  achieve reliable ECG fragments and to divide them beat by beat. Then, for each heartbeat, 3D FMM<sub>ecg</sub> provides FMM parameter estimates. For each wave, the series of parameter values, corresponding to consecutive beats, which can be summarized to get average patterns as well as the changes in the patterns over time.

### Data Preprocessing

3D FMM<sub>ecg</sub> code incorporates a standard preprocessing for 1-12 ECG leads including baseline correction, QRS detection and ECG segmentation.

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

PTB-XL is a large dataset of 21837 clinical 12-lead ECGs of 10 second length annotated by two cardiologists with diagnostic labels, based on SCP-ECG statements, and the likelihood information for the  statements [4]. PTB-XL database has been analysed using 3D FMM<sub>ecg</sub> in [1] for patients with likelihood $\geq 80$. In particular, we analysed 9055 patients from PTB-XL labelled as NORM, i.e. with non-pathological ECGs. 

The normal percentile ranges ( $5th$, $95th$ ) of FMM parameters, and indices defined from them, calculated from NORM patients are specially useful for identifyng noisy and/or pathological ECG patterns, as those subjects with values out such ranges. The FMM-based indices for which these ranges were computed are: 

* $Var_J$: measure of the relative relevance of wave $J$, for $J= P, Q, R, S, T$. For a given lead, it is calculated as the variability the wave $J$ explains, see [5] for details.
* $R^2:$ measure of the accuracy of the model across leads, see [5] for details.
* RR: duration of R-R interval from QRS annotations in milisecons (ms).
* disPQ, disQS, and disQT: difference in ms between $\alpha_P$ and $\alpha_Q;$, $\alpha_Q$ and $\alpha_S$; and $\alpha_S$ and $\alpha_T$, respectively.

Normal percentile ranges of FMM parameters and these indices across NORM patients in PTB-XL are given in Tables 1-5 for the medians (Me) and in Tables 6-10 for the coefficient of variation (Cv). Cv for angular parameters was defined in the Suplementary Material of [5]. 


||$5th_P$|$95th_P$|$Mean_P$|$Sd_P$|$5th_Q$|$95th_Q$|$Mean_Q$|$Sd_Q$|$5th_R$|$95th_R$|$Mean_R$|$Sd_R$|$5th_S$|$95th_S$|$Mean_S$|$Sd_S$|$5th_T$|$95th_T$|$Mean_T$|$Sd_T$|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---
I_A|22.002|73.627|44.938|15.852|57.324|262.948|138.689|64.943|154.025|638.104|371.694|149.128|34.932|218.233|107.977|57.965|51.812|192.935|113.685|43.848
II_A|30.253|109.85|66.787|24.282|46.149|254.907|132.262|66.475|183.181|722.084|420.162|166.09|42.612|254.059|128.446|66.519|63.437|220.038|130.779|48.231
V1_A|13.868|59.663|33.515|14.337|61.708|300.483|159.988|76.569|173.454|757.599|435.613|181.521|67.354|372.525|196.125|93.977|19.53|134.936|63.486|38.089
V2_A|15.96|72.734|39.336|17.644|103.525|527.64|273.484|133.914|217.064|1189.753|616.442|303.316|147.713|636.997|362.287|149.064|60.272|450.34|228.606|119.734
V3_A|17.814|70.079|39.857|16.37|105.404|479.718|258.415|118.248|203.131|1061.654|573.514|270.119|166.47|644.52|373.345|149.792|84.828|464.813|247.689|117.127
V4_A|19.753|67.329|40.85|15.005|106.309|508.348|263.392|127.063|336.559|1230.165|726.805|277.21|129.452|584.903|318.72|144.369|82.262|398.297|211.929|98.907
V5_A|20.243|64|39.969|13.754|102.544|448.654|236.717|110.61|350.096|1115.739|692.33|232.764|73.632|423.253|220.787|110.034|74.012|303.489|170.084|73.399
V6_A|19.225|60.177|37.595|12.606|75.297|333.714|178.898|83.4|272.19|873.607|544.397|184.67|46.298|297.477|149.211|80.178|56.324|229.983|129.737|55.291
I_Beta|3.065|4.617|3.846|0.475|0.082|2.246|1.237|0.64|2.977|4.581|3.525|0.49|3.734|0.823|5.426|1.042|2.744|3.544|3.215|0.252
II_Beta|3.378|4.808|4.069|0.463|6.267|2.327|1.23|0.697|2.935|4.467|3.455|0.47|3.493|0.889|5.454|1.085|2.655|3.47|3.122|0.263
V1_Beta|4.854|1.37|0.115|0.812|1.858|4.816|3.31|0.91|5.413|0.344|6.078|0.389|5.976|2.511|0.666|0.748|3.583|0.255|5.06|1.11
V2_Beta|3.802|0.537|5.453|0.917|1.381|4.153|2.559|0.866|3.922|0.106|5.518|0.752|5.723|1.822|0.41|0.72|3.243|4.009|3.608|0.315
V3_Beta|3.598|6.184|4.906|0.79|0.96|3.259|1.944|0.704|3.195|6.066|4.415|0.927|5.134|1.158|0.06|0.728|3.065|3.756|3.432|0.246
V4_Beta|3.451|5.574|4.504|0.657|0.397|2.469|1.443|0.613|3.05|5.197|3.755|0.643|5.015|0.981|6.144|0.702|2.916|3.607|3.29|0.229
V5_Beta|3.316|5.161|4.235|0.567|0.075|2.178|1.209|0.619|3.026|4.67|3.56|0.503|4.709|0.884|5.882|0.757|2.789|3.524|3.195|0.236
V6_Beta|3.244|4.958|4.071|0.528|6.165|2.069|1.088|0.647|2.956|4.295|3.412|0.416|3.931|0.607|5.359|0.887|2.659|3.471|3.126|0.259
ALL_Alpha|4.417|5.152|4.847|0.224|5.506|5.637|5.609|0.041|5.745|5.876|5.765|0.039|5.862|6.005|5.889|0.048|0.954|1.766|1.336|0.254
ALL_Omega|0.067|0.219|0.135|0.047|0.024|0.064|0.039|0.016|0.026|0.044|0.035|0.008|0.019|0.051|0.031|0.01|0.124|0.288|0.181|0.052

Table 1: $5th$, $95th$, Mean and standard deviation (Sd) for the Me of FMM parameters. Subindices in the header denote waves. ALL is used to denote a common estimation across leads. 

||$5th_P$|$95th_P$|$Mean_P$|$Sd_P$|$5th_Q$|$95th_Q$|$Mean_Q$|$Sd_Q$|$5th_R$|$95th_R$|$Mean_R$|$Sd_R$|$5th_S$|$95th_S$|$Mean_S$|$Sd_S$|$5th_T$|$95th_T$|$Mean_T$|$Sd_T$|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---
I|0.006|0.068|0.028|0.022|0.014|0.183|0.071|0.057|0.211|0.822|0.545|0.19|0.004|0.094|0.034|0.03|0.05|0.548|0.254|0.154
II|0.006|0.162|0.056|0.054|0.009|0.135|0.054|0.042|0.2|0.801|0.53|0.186|0.005|0.133|0.047|0.044|0.06|0.519|0.254|0.142
V1|0.001|0.055|0.017|0.021|0.009|0.186|0.069|0.06|0.41|0.892|0.72|0.152|0.013|0.207|0.075|0.071|0.005|0.243|0.072|0.082
V2|0.001|0.026|0.008|0.01|0.01|0.148|0.062|0.047|0.131|0.806|0.488|0.208|0.015|0.222|0.09|0.072|0.045|0.671|0.322|0.192
V3|0.001|0.034|0.01|0.012|0.011|0.136|0.06|0.041|0.096|0.698|0.394|0.182|0.021|0.247|0.113|0.073|0.102|0.721|0.388|0.186
V4|0.001|0.033|0.01|0.012|0.013|0.18|0.07|0.053|0.263|0.772|0.532|0.156|0.018|0.185|0.081|0.054|0.076|0.53|0.272|0.139
V5|0.001|0.032|0.01|0.012|0.016|0.213|0.076|0.064|0.365|0.833|0.625|0.144|0.007|0.114|0.046|0.036|0.057|0.404|0.209|0.109
V6|0.002|0.041|0.013|0.015|0.014|0.201|0.07|0.063|0.396|0.865|0.669|0.144|0.004|0.082|0.031|0.026|0.042|0.374|0.181|0.104
ALL|0.001|0.068|0.019|0.028|0.012|0.173|0.066|0.055|0.199|0.843|0.563|0.198|0.007|0.182|0.065|0.061|0.024|0.572|0.244|0.168

Table 2: $5th$, $95th$, Mean and standard deviation (Sd) for the Me of $Var_J$. Subindices in the header denote waves. ALL is used to denote a common estimation across leads. 

||$5th$|$95th$|$Mean$|$Sd$|
|---|---|---|---|---
I|0.875|0.984|0.949|0.047
II|0.901|0.986|0.958|0.033
V1|0.921|0.987|0.965|0.026
V2|0.957|0.989|0.977|0.012
V3|0.955|0.989|0.977|0.012
V4|0.963|0.991|0.981|0.01
V5|0.963|0.992|0.982|0.01
V6|0.958|0.991|0.98|0.016
ALL|0.93|0.99|0.971|0.027

Table 3: $5th$, $95th$, Mean and standard deviation (Sd) for the Me of $R^2$. ALL is used to denote a common estimation across leads. 

||$5th$|$95th$|$Mean$|$Sd$|
|---|---|---|---|---
RR (ms)|636|1152|881.487|157.947

Table 4: $5th$, $95th$, Mean and standard deviation (Sd) for the Me of $RR$.

||$5th$|$95th$|$Mean$|$Sd$|
|---|---|---|---|---
disPQ (ms)|72|142|104|25
disQS (ms)|28|49|39|7
disQT (ms)|234|322|277|27

Table 5: $5th$, $95th$, Mean and standard deviation (Sd) for the Me of disPQ, disQS and disQT.

||$5th_P$|$95th_P$|$Mean_P$|$Sd_P$|$5th_Q$|$95th_Q$|$Mean_Q$|$Sd_Q$|$5th_R$|$95th_R$|$Mean_R$|$Sd_R$|$5th_S$|$95th_S$|$Mean_S$|$Sd_S$|$5th_T$|$95th_T$|$Mean_T$|$Sd_T$|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---
I_A|10.928|50.222|25.679|12.944|16.625|65.51|38.566|15.301|5.178|28.473|13.622|7.901|18.704|88.235|47.692|21.649|3.878|20.446|9.745|6.509
II_A|6.572|42.287|18.954|11.837|18.086|75.997|43.618|17.979|4.193|27.188|12.244|7.996|17.3|86.689|45.699|21.449|3.12|18.305|8.352|5.879
V1_A|12.572|67.382|35.498|17.367|12.388|61.078|31.909|15.415|5.688|42.581|19.05|12.088|10.021|63.212|30.458|17.009|6.372|47.676|20.457|13.295
V2_A|17.189|73.466|41.816|17.981|8.888|49.06|24.412|13.328|6.309|42.357|20.512|11.865|8.985|53.119|27.074|14.267|3.094|24.072|9.802|8.266
V3_A|15.474|72.445|40|18.585|9.52|50.758|25.921|13.352|5.698|44.224|19.446|12.789|10.311|57.138|29.84|14.897|2.535|17.316|7.34|5.97
V4_A|12.735|65.691|34.808|17.259|13.726|61.661|35.786|14.938|4.749|28.179|13.421|8.021|13.21|61.258|33.624|15.293|2.702|15.836|7.192|5.187
V5_A|12.317|59.843|31.538|15.9|16.806|63.696|39.196|14.605|3.964|25.109|11.382|7.146|17.112|81.028|43.112|20.049|2.85|16.808|7.593|5.841
V6_A|11.459|59.432|30.312|15.642|17.752|67.63|41.044|15.318|3.815|24.12|10.805|7.001|19.68|90.592|49.819|21.989|3.187|20.392|8.794|6.509
I_Beta|24.282|56.38|40.628|9.816|22.802|53.785|38.539|9.376|27.275|55.85|41.234|8.676|19.989|67.366|41.528|14.304|7.332|40.707|24.123|10.863
II_Beta|19.621|54.962|36.711|10.757|21.682|56.611|39.361|10.586|25.521|59.219|42.017|10.291|19.705|68.487|41.441|14.651|7.537|42.071|23.81|10.621
V1_Beta|36.747|77.158|58.934|12.287|36.378|74.292|57.668|11.536|44.184|81.63|65.828|11.398|20.237|69.28|45.572|14.684|12.383|88.347|52.855|26.236
V2_Beta|22.537|64.855|43.456|12.919|21.326|66.306|43.974|13.897|27.265|75.437|51.834|15.173|18.414|61.936|39.665|13.168|7.632|39.009|22.516|10.781
V3_Beta|19.607|56.621|36.697|11.286|17.981|54.591|34.163|10.97|22.667|60.293|37.997|11.3|16.987|54.994|34.5|11.658|6.154|36.465|20.641|9.944
V4_Beta|18.848|53.499|35.003|10.637|17.987|50.602|33.753|9.903|21.81|54.682|37.175|9.859|16.43|52.005|32.634|10.812|5.845|40.212|21.338|10.951
V5_Beta|19.591|54.976|36.624|10.683|20.226|54.156|37.45|10.38|23.669|56.971|39.993|10.193|17.146|52.715|33.541|10.871|6.673|40.579|23.225|11.083
V6_Beta|19.871|52.473|35.739|9.762|21.267|53.164|37.794|9.689|25.248|53.53|39.259|8.718|19.475|61.416|38.247|12.763|7.265|36.038|22.317|9.418
ALL_Alpha|1.47|12.908|4.92|3.623|0.086|4.157|1.878|1.276|0.094|3.884|1.647|1.213|0.078|4.211|1.852|1.384|0.936|6.667|3.013|2.041
ALL_Omega|11.563|57.799|30.612|14.601|14.499|84.667|40.736|22.449|2.034|21.064|12.45|6.793|10.663|59.77|31.103|17.876|1.874|20.523|9.403|6.721

Table 6: $5th$, $95th$, Mean and standard deviation (Sd) for the Cv of FMM parameters. Subindices in the header denote waves. ALL is used to denote a common estimation across leads. 

||$5th_P$|$95th_P$|$Mean_P$|$Sd_P$|$5th_Q$|$95th_Q$|$Mean_Q$|$Sd_Q$|$5th_R$|$95th_R$|$Mean_R$|$Sd_R$|$5th_S$|$95th_S$|$Mean_S$|$Sd_S$|$5th_T$|$95th_T$|$Mean_T$|$Sd_T$|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---
I|24.47|103.221|55.216|25.835|28.63|114.317|63.618|27.06|3.288|35.141|14.643|10.624|38.751|149.583|85.326|34.671|7.146|45.511|20.862|13.807
II|16.146|89.965|43.251|24.441|30.183|129.517|71.119|31.468|3.086|35.655|13.908|12.008|35.113|146.066|81.823|34.87|6.161|40.8|18.411|12.6
V1|30.711|140.92|76.444|34.462|25.158|115.452|61.67|28.747|2.296|30.834|11.06|12.429|21.059|120.989|61.348|31.819|13.429|112.967|49.096|32.214
V2|40.24|163.879|91.638|38.417|20.713|106.09|54.982|27.48|3.505|43.296|17.507|14.705|19.522|109.739|55.8|29.666|4.168|48.228|18.799|17.564
V3|39.43|171.504|94.107|41.527|24.406|117.053|63.021|29.162|4.628|51.949|20.722|16.796|20.976|105.943|55.482|28.336|3.582|37.546|14.807|12.983
V4|36.778|160.112|86.587|38.683|26.432|115.994|64.826|28.501|3.274|27.893|12.32|9.024|27.118|113.591|63.012|27.332|5.098|37.275|16.478|12.14
V5|34.341|145.713|78.247|35.316|27.923|113.791|64.821|26.793|2.705|26.212|11.412|7.999|35.646|141.952|79.514|33.447|6.195|40.579|18.453|12.084
V6|30.973|132.269|71.039|32.86|29.613|118.508|67.609|28.001|2.423|27.487|11.196|8.643|40.942|155.563|89.903|35.816|7.465|47.822|21.729|14.257

Table 7: $5th$, $95th$, Mean and standard deviation (Sd) for the Cv of $Var_J$. Subindices in the header denote waves. 

||$5th$|$95th$|$Mean$|$Sd$|
|---|---|---|---|---
I|0.396|4.436|1.714|1.795
II|0.367|3.703|1.514|1.385
V1|0.332|4.039|1.585|2.103
V2|0.264|2.478|1.056|0.884
V3|0.255|2.565|1.096|1.064
V4|0.208|2.17|0.923|0.904
V5|0.196|2.251|0.913|0.998
V6|0.217|2.547|1.033|1.63

Table 8: $5th$, $95th$, Mean and standard deviation (Sd) for the Cv of $R^2$.  

||$5th$|$95th$|$Mean$|$Sd$|
|---|---|---|---|---
RR (ms)|0.766|7.788|3.168|2.341

Table 9: $5th$, $95th$, Mean and standard deviation (Sd) for the Cv of $RR$.

||$5th$|$95th$|$Mean$|$Sd$|
|---|---|---|---|---
disPQ (ms)|5|52|19|15
disQS (ms)|2|36|17|12
disQT (ms)|1|6|3|2

Table 10: $5th$, $95th$, Mean and standard deviation (Sd) for the Cv of disPQ, disQS and disQT.


## References

[1] Rueda, C., Rodríguez-Collado, A., Fernández, I., Canedo, C., Ugarte, M. D., & Larriba, Y. (2022). A unique cardiac electrocardiographic 3D model. Toward interpretable AI diagnosis. iScience, 25(12), 105617. https://doi.org/10.1016/j.isci.2022.105617

[2] Fernández, I., Rodríguez-Collado, A., Larriba, Y., Lamela, A., Canedo, C., and Rueda, C. (2021). FMM: rhythmic patterns modeling by FMM models. R package version 0.3.0.

[3] Fernández, I., Rodríguez-Collado, A., Larriba, Y., Lamela, A., Canedo, C., and Rueda, C. (2022). FMM: an R package for modeling rhythmic patterns in oscillatory systems. R Journal. 14, 361–380.

[4] Wagner, P., Strodthoff, N., Bousseljot, RD. et al. (2020) PTB-XL, a large publicly available electrocardiography dataset. Scientific Data, 7, 154. https://doi.org/10.1038/s41597-020-0495-6

[5] Rueda, C.,  Larriba, Y., & Lamela, A. (2021). The hidden waves in the ECG uncovered revealing a sound automated interpretation method. Scientific Reports, 11, 3724. https://doi.org/10.1038/s41598-021-82520-w

