#!/bin/bash
## IDR thresholds default: 0.05 (similar to idea of FDR)

SAMPLE=pcrFree.st105_foxh1  ####change this everytime for new samples####

## Nt= rep1 vs rep2 true reps compare
idr --samples sorted.${SAMPLE}_rep1_peak_peaks.narrowPeak sorted.${SAMPLE}_rep2_peak_peaks.narrowPeak --input-file-type narrowPeak --rank signal.value --output-file ${SAMPLE}_rep1_vs_rep2_idr --plot --log-output-file ${SAMPLE}_rep1_vs_rep2_idr_log
##filter out peaks idr < 0.05
awk '{if($5 >= 540) print $0}' ${SAMPLE}_rep1_vs_rep2_idr > ${SAMPLE}_rep1_vs_rep2_idr_0.05_filtered.narrowPeak

## Np= pooled rep1 and rep2 then split for 2 pseudo reps to compare
idr --samples sorted.merged_${SAMPLE}_00_peak_peaks.narrowPeak sorted.merged_${SAMPLE}_01_peak_peaks.narrowPeak --input-file-type narrowPeak --rank signal.value --output-file ${SAMPLE}_pooled_pseudo00_vs_pseudo01_idr --plot --log-output-file ${SAMPLE}_pooled_pseudo00_vs_pseudo01_idr_log
##filter out peaks idr < 0.05
awk '{if($5 >= 540) print $0}' ${SAMPLE}_pooled_pseudo00_vs_pseudo01_idr > ${SAMPLE}_pooled_pseudo00_vs_pseudo01_idr_0.05_filtered.narrowPeak

## N1= rep1 splits into 2 pseudo samples to compare
idr --samples sorted.${SAMPLE}_rep1_00_peak_peaks.narrowPeak sorted.${SAMPLE}_rep1_01_peak_peaks.narrowPeak --input-file-type narrowPeak --rank signal.value --output-file ${SAMPLE}_rep1_pseudo00_vs_pseudo01_idr --plot --log-output-file ${SAMPLE}_rep1_pseudo00_vs_pseudo01_idr_log
##filter out peaks idr < 0.05
awk '{if($5 >= 540) print $0}' ${SAMPLE}_rep1_pseudo00_vs_pseudo01_idr > ${SAMPLE}_rep1_pseudo00_vs_pseudo01_idr_0.05_filtered.narrowPeak

## N2= rep2 splits into 2 pseudo samples to compare
idr --samples sorted.${SAMPLE}_rep2_00_peak_peaks.narrowPeak sorted.${SAMPLE}_rep2_01_peak_peaks.narrowPeak --input-file-type narrowPeak --rank signal.value --output-file ${SAMPLE}_rep2_pseudo00_vs_pseudo01_idr --plot --log-output-file ${SAMPLE}_rep2_pseudo00_vs_pseudo01_idr_log
##filter out peaks idr < 0.05
awk '{if($5 >= 540) print $0}' ${SAMPLE}_rep2_pseudo00_vs_pseudo01_idr > ${SAMPLE}_rep2_pseudo00_vs_pseudo01_idr_0.05_filtered.narrowPeak
#### quality control for reps
####(N1/N2 >= 2)& (Np/Nt >= 2) then it is low quality rep or low reproducibility
