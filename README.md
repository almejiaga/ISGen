# ISGen
[![Build Status](https://travis-ci.org/DomNelson/ISGen.svg?branch=master)](https://travis-ci.org/DomNelson/ISGen)

Importance sampling in large genealogies


#   1) Prepare your input files

For running ISgen you will need a few files:
## 1.1) ped file

The ped file is a standard plink ped file, with at least three columns:
```
ind father mother sex
1 11 12 1
2 15 14 2
3 15 14 2
11 102 101 1
12 0 0 2
13 102 101 1
14 0 0 2
15 103 104 1
16 103 104 2
18 105 106 2
19 105 106 2
20 107 108 2
21 107 108 1
101 0 0 2
102 202 201 1
103 0 0 1
104 202 201 2
105 202 201 1
106 0 0 2
107 202 201 1
108 0 0 2
201 0 0 2
202 0 0 1
```
## 1.2) List of heterozygous and/or homozygous IDs

The file is a one column file with the IDs that correspond to Het IDs, like this:
```
1
2
3
```
If you also have homozygous IDs, they should be listed in another file

With this two files you can run the climbing part of this software. If you also want to estimate the regional frequencies, you need two additional files:

## 1.3) map folder

Note: this script recognizes vcfs files that end in .vcf.gz

---
#   2) Generate your genetic map file in increasing genetic distance for ARG_needle
for this section, we will use the .bim file created in the last step. the R script is designed to work automatically but we need to create two folders in order to the script to work properly

## 2.1) map folder
We have to create a map folder that will store the .bim files for each chromosome

## 2.2) plink_maps folder
We have to create a plink_maps folder that stores the genetic map files for each chromosome in plink format. This files can be found in different resources, such as: put the link here. 

## 2.3) Generate a list of excluded variants from the analysis
Now that we have our files in the respective folders, we just have to run the script in this way in linux:

Rscript geneticmapfilepreparation.R

After this, you will get a list of variants that are going to be excluded from the analysis only for the fact that the genetic distance between them is to small that ARG_needle is not able to say that you have the variants in an increasing number. 

# 3) excluding the variants from your vcf to be able to run ARG_needle
Run the script:

bash filteringagain.sh

With this script you will create a folder named "ARG_files" and the .haps and .samples files will be stored there.
