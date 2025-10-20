# Bothrops insularis genome
<!--- [![Published in GBE](https://img.shields.io/badge/published%20in-GBE-blue)](https://doi.org/10.1093/molbev/msaf058) --->

[![Data available in the Fisghare](https://img.shields.io/badge/data%20available%20in%20the-figshare-red)](https://figshare.com/projects/Bothrops_insularis_genome/237995)

<!---
This repository contains commands and scripts used in the manuscript "The golden lancehead genome reveals distinct selective processes acting on venom genes of an island endemic snake" published in *Genome Biology and Evolution*.

All datasets used in the present study are detailed in the Supplementary file of the published manuscript.
--->

The datasets are available in the NCBI under the project number [PRJNA679826](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA679826/).

## Genome assembly
The genome assembly pipeline is described in this turorial: "[Tutorial to perform chromosome-level genome assembly using HiFi and HiC data](https://github.com/pedronachtigall/HI-genome-assembly-pipeline)".

## Genome annotation
### Repeat annotation
The repeat annotation was performed using [RepeatModeler2](https://github.com/Dfam-consortium/RepeatModeler), to generate a *de novo* species-specific library, and [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker), to perform the annotation. We complement the species-specific TE library using a curated TE library of snakes previosuly published <sup>[Castoe et al., 2013](https://doi.org/10.1073/pnas.1314475110)</sup>.

The pipeline with commands and scripts used to perform the repeat annotation is decribed in the following tutorial: https://github.com/pedronachtigall/Repeat-annotation-pipeline

For this step, we used the primary genome assembly to perform the repeat annotation and the soft-masked primary genome assembly as the source for gene annotation.

### Gene annotation
The gene annotation was performed using [GALBA](https://github.com/Gaius-Augustus/GALBA). We followed the commands decribed in the following tutorial: "[Gene-annotation-pipeline](https://github.com/pedronachtigall/Gene-annotation-pipeline)".

For this step, we used the soft-maked primary assembly to perform the gene annotation using GALBA with the integration of RNA-seq data derived from several tissues to be used as transcript evidence.

The mitochondrial genome was annotated using [MitoZ](https://github.com/linzhi2013/MitoZ).
```
MitoZ.py annotate --genetic_code auto --clade Chordata --outprefix mitogenome_annotation --thread_number 20 --fastafile Binsularis_mitogenome.fasta
```

### Toxin annotation
We used [ToxCodAn-Genome](https://github.com/pedronachtigall/ToxCodAn-Genome) to annotate toxins.

For this step, we used the primary assembly, the venom-gland transcriptome data annotated using [ToxCodAn](https://github.com/pedronachtigall/ToxCodAn) and the [Viperidae](https://raw.githubusercontent.com/pedronachtigall/ToxCodAn-Genome/main/Databases/Viperidae_db.fasta) database. We used the venom-gland transcriptome data of the individual used to seqeunce the genome (SRR32358140).

First, we assembled the venom-gland transcriptome. For the transcriptome assembly, we used the "[TRassembly.py](https://raw.githubusercontent.com/pedronachtigall/ToxCodAn-Genome/refs/heads/main/bin/TRassembly.py)" script from ToxCodAn-Genome, which uses the genome-guided mode of Trinity and StringTie and the *de novo* mode of Trinity and rnaSPAdes.
```
#trim reads
trim_galore --paired --phred33 --length 75 -q 25 --stringency 1 -e 0.1 -o SRR32358140_tg SRR32358140_R1.fastq.gz SRR32358140_R2.fastq.gz

#assembly transcripts
TRassembly.py -g Cadam_primary_chromosomes.fasta -r SRR32358140_R1_val_1.fastq.gz,SRR32358140_R2_val_2.fastq.gz -c 20 -M 20G --output SRR32358140_TRassembly
```

Annotate the transcriptome using ToxCodAn.
```
wget https://github.com/pedronachtigall/ToxCodAn/raw/refs/heads/master/models.zip
unzip models.zip
toxcodan.py -c 10 -t SRR32358140_TRassembly/transcripts.fasta -m models -o SRR32358140_ToxCodAn -s SRR32358140

cat SRR32358140_ToxCodAn/SRR32358140_Toxins_cds_RedundancyFiltered.fasta SRR32358140_ToxCodAn/SRR32358140_PutativeToxins_cds_SPfiltered.fasta > Binsularis_VG_toxins.toxcodan.fasta
```

Then, run ToxCodAn-Genome setting the annotated toxins from the venom-gland transcriptome and the Viperidae database.
```
wget https://raw.githubusercontent.com/pedronachtigall/ToxCodAn-Genome/main/Databases/Viperidae_db.fasta

toxcodan-genome.py -c 10 -g Cadam_primary_chromosomes.fasta -C Binsularis_VG_toxins.toxcodan.fasta -d Viperidae_db.fasta
```

The toxin genes were manually reviewed and curated following the "[Checking annotations](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Guide#checking-annotations)" section of the ToxCodAn-Genome's guide to ensure a comprehensive toxin annotation.

Additionally, we compared our toxin-annotated sequences to the ToxProt database, which comprises a curated and reviewed set of toxin sequences, and to the expressed sequence tag (EST) of the venom gland of *B. insularis* (NCBI accessions from BM401391 to BM402067) using BLAST.
```
ADD CODE
```

## Estimate expression level
We mapped the transcriptomic reads of several tissues against the *B. insularis* genome using [STAR](https://github.com/alexdobin/STAR/releases) and retrieved the counts using [featureCounts](https://subread.sourceforge.net/featureCounts.html).

The raw data is listed below:
| Sample ID | Data type | NCBI accession |
| :-------- | :-------: | :------------: | 
| SB1851_VG_rna  | RNA-seq | SRR32358140 |
| SB1851_HG_rna  | RNA-seq | SRR32358139 |
| SB1851_ILG_rna  | RNA-seq | SRR32358138 |
| SB1851_Pancreas_rna  | RNA-seq | SRR32358137 |
| SB1851_Muscle_rna  | RNA-seq | SRR32358136 |
| SB1851_Heart_rna  | RNA-seq | SRR32358151 |
| SB1851_Brain_rna  | RNA-seq | SRR32358150 |
| SB1851_Spleen_rna  | RNA-seq | SRR32358149 |
| SB1851_Kidney_rna  | RNA-seq | SRR32358148 |
| SB1851_Ovary_rna  | RNA-seq | SRR32358147 |
| SB0825_VG_rna  | RNA-seq | SRR32358146 |
| SB1543_VG_rna  | RNA-seq | SRR32358145 |
| SB1852_VG_rna  | RNA-seq | SRR32358144 |

```
ADD CODE
```

## Analysis of genetic variation in toxin genes
We used whole genome sequencing data from eight specimens of *B. insularis* to check for the genetic variation of toxin genes.

The raw data is listed below:
| Sample ID | Data type | NCBI accession |
| :-------- | :-------: | :------------: | 
| SB1851    | WGS       | SRR32358141    |
| SB0825    | WGS       | SRR32361269    |
| SB1852    | WGS       | SRR32361268    |
| SB1853    | WGS       | SRR32361267    |
| SB1854    | WGS       | SRR32361266    |
| SB1855    | WGS       | SRR32361265    |
| SB1856    | WGS       | SRR32361264    |
| SB1858    | WGS       | SRR32361263    |

```
ADD CODE
```

### Demographic inference
We used [SMC++](https://github.com/popgenmethods/smcpp) to estimate effective population size.
```
ADD CODE
```

## Comparison to *Bothrops jararaca*
We compared the toxin genes from *B. insularis* to the toxin anntoation from *B. jararaca* genome ([Almeida et al., 2021](https://doi.org/10.1073/pnas.2015159118); GCA_018340635.1).

### Orthology inference
We used Orthofinder to infer orthology of toxin genes.
```
ADD CODE
```

### Gene tree
When needed (e.g., to analyze the multi-loci toxin genes, such as SVMPs, SVSPs, CTLs, and PLA2s), we inferred the gene tree using their CDSs. We used [MAFFT](https://github.com/GSLBiotech/mafft) and [IQ-TREE](https://github.com/Cibiv/IQ-TREE) for alignment of sequences and infer gene tree, respectively. The input files, which consist of a fasta file with sequences from the target gene family, were set manually.
```
mafft --auto GENE.fasta > GENE.ALIGNED.fasta
iqtree -s GENE.ALIGNED.fasta -m TEST -bb 1000 -alrt 1000
```

The trees were inspected and adjusted using [FigTree](https://github.com/rambaut/figtree/).

## Macrosyntenic analysis
We compared the *B. insularis* chromosome genome assembly with other viper species with a near chromosome-level genome assembly available in the literature (i.e., *Crotalus adamanteus*, *Crotalus viridis*, *Sistrurus catenatus*, *Gloydius shedaoensis*, *Cerastes gasperettii*, and *Vipera ursinii*). To do this, we used [GENESPACE](https://github.com/jtlovell/GENESPACE) and followed its guide to compare genomes at gene-level.

## Cite
If you follow the pipeline and/or use any of the scripts in this repository, please cite:

Nachtigall et al. (in prep)

```
```
