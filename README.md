# Mutation Calling Pipeline

This is the pipeline for DNA mutation calling. 

# Outline

<!-- TOC -->

- [Mutation Calling Pipeline](#mutation-calling-pipeline)
- [Outline](#outline)
- [Introduction](#introduction)
- [Copyright](#copyright)
- [Running Environment](#running-environment)
- [Useage](#useage)
    - [0. Introduction](#0-introduction)
    - [1. Data Downloading and Quality Control](#1-data-downloading-and-quality-control)
    - [2. Data Cleaning and its Quality Control, Alignment and Recabliration, Mutation Calling, Mutation Filtering](#2-data-cleaning-and-its-quality-control-alignment-and-recabliration-mutation-calling-mutation-filtering)
- [Reference](#reference)

<!-- /TOC -->

# Introduction

This pipeline can be used for analyzing the next generation sequencing data. It includes six steps:
1. Data Downloading (*This step can be skipped if the sequencing data is already on hand.*)

2. Quality Control

3. Data Cleaning and its Quality Control

4. Alignment and Recabliration

5. Mutation Calling

6. Mutation Filtering

The corresponding running scripts are provided. For their useage, please refer to the later content.

# Copyright

Copyright (C) 2020 Edison Jianning KANG

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

![ [Large GPLv3 logo with “Free as in Freedom”] ](https://www.gnu.org/graphics/gplv3-with-text-136x68.png)

# Running Environment

1. FastQC >=0.11.9

2. CutAdapt >=2.10

3. BWA >=0.7.17-r1188

4. GATK, both v3 and v4

   The output of `gatk --version`:

   ```
   The Genome Analysis Toolkit (GATK) v4.1.4.1
   HTSJDK Version: 2.21.0
   Picard Version: 2.21.2
   ```

   GATK 3 >=`3.8-1-0-gf15c1c3ef`

6. MuTect1 >=3.1-0-g72492bb

7. VarScan >=2.4.2

8. Manta >=1.6.0.centos6_x86_64

9. Strelka >=2.9.10

10. IBM Aspera

    ```
    Aspera Connect version 3.9.8.176272
    ascp version 3.9.1.168302
    ```
    
11. ANNOVAR

    ```
    Version: $Date: 2019-10-24 00:05:27 -0400 (Thu, 24 Oct 2019) $
    ```
    
12. R>=4.0.2

13. VCFTools>=0.1.17

14. CrossMap>=0.4.1

15. gzip>=1.10

16. Unzip>=6.00

# Useage
## 0. Introduction
All the scripts and an example directory structure for running those scripts have been added to this repository. To use the pipeline, you can directly fork or clone the repository to the local computer. 

(**Notice: Any further change may be added to the pipeline. If you want to keep using the latest one, it is better to fork the repository rather than clone it.**)

The structure of this repository is illustrated below:
```
    - codes
    It stores all the scripts and an example directory structure
        - 1_ORIGIN
        - 2_FASTQC
        - 3_CUTADAPT
        - 4_BWA_GATK
        - 5_VARIANT_CALLING
        - 6_FILTERING
        - log

        Above directories are an example directory strcture. Running results of each steps will be stored in corresponding dicretory. For log directory, it will store the log file of all the steps.

          1_Aspera_FASTQC.sh
          All.sh

        Above are all the scripts needed for running.

          2_PLOT_IN_BRIEF.R
          3_1_FOR_REFERENCE.R
          3_2_DEPENDT_FILE_add_fs_annotation.R
          add_readct_shimmer.pl
          BEFORE_ANNOVAR_FILTERING.R
          AFTER_ANNOVAR_FILTERING.R

        Above are all the dependent files that are needed for running the analysis.

          hg382hg19.sh
          
        Above is the scripts that used to convert mutations from hg38 to hg19.

    - LICENSE
    It is the lincense file for this pipeline

    - README.md
    It is the ReadMe file that shown here.
```

## 1. Data Downloading and Quality Control
To perform the first two steps, you can run the script [`1_Aspera&FASTQC.sh`](https://github.com/Edison19991109/Mutation_Calling/blob/master/1_Aspera%26FastQC.sh). Downloading is conducted by using the `IBM Aspera` and an example downloading ftp link is provided. Quality control is conducted by using `FastQC`. To run the script, sample header and sample ID are needed.

*(Detail version information of all the programs that are mentioned here are provided in an early section.)*

## 2. Data Cleaning and its Quality Control, Alignment and Recabliration, Mutation Calling, Mutation Filtering
To run the last four steps, you can run the script [`All.sh`](https://github.com/Edison19991109/Mutation_Calling/blob/master/All.sh). Data cleaning is done through `Cutadapt` and an example command of running Cutadapt is provided. Qulity control is also conducted by using `FastQC` (ommited in this version of the script). Alignment is conducted by `BWA-MEM` and recabliration is conducted by `GATK Best Practice` pipeline. Mutation calling is conducted parallelly by 6 different mutation calling programs and all the results are combined together and then filtered by using R scripts and the `Annovar` annotation result.


*(Detail version information of all the programs that are mentioned here are provided in an early section.)*

# Reference
https://github.com/tianshilu/QBRC-Somatic-Pipeline
