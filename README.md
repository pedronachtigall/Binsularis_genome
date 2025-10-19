# Bothrops insularis genome
[![Published in GBE](https://img.shields.io/badge/published%20in-GBE-blue)](https://doi.org/10.1093/molbev/msaf058)
[![Data available in the Fisghare](https://img.shields.io/badge/data%20available%20in%20the-figshare-red)](https://figshare.com/projects/Bothrops_insularis_genome/237995)

This repository contains commands and scripts used in the manuscript "The golden lancehead genome reveals distinct selective processes acting on venom genes of an island endemic snake" published in *Genome Biology and Evolution*.

All datasets used in the present study are detailed in the Supplementary file of the published manuscript.

## Genome assembly
The genome assembly pipeline is described in this turorial: "[Tutorial to perform chromosome-level genome assembly using HiFi and HiC data](https://github.com/pedronachtigall/HI-genome-assembly-pipeline)".

## Genome annotation
### Repeat annotation
The repeat annotation was performed using [RepeatModeler2](https://github.com/Dfam-consortium/RepeatModeler), to generate a *de novo* species-specific library, and [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker), to perform the annotation. We complement the species-specific TE library using a curated TE library of snakes previosuly published <sup>[Castoe et al., 2013](https://doi.org/10.1073/pnas.1314475110)</sup>.

The pipeline with commands and scripts used to perform the repeat annotation is decribed in the following tutorial: https://github.com/pedronachtigall/Repeat-annotation-pipeline

For this step, we used the primary genome assembly to perform the repeat annotation and the soft-masked primary genome assembly as the source for gene annotation.

### Gene annotation
