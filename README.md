# ISGen
[![Build Status](https://travis-ci.org/DomNelson/ISGen.svg?branch=master)](https://travis-ci.org/DomNelson/ISGen)

Importance sampling in large genealogies


#   1) Prepare your input files

For running ISgen you will need a four files:

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

With these two files you can run the climbing part of this software. If you also want to estimate the regional frequencies, you need two additional files:

## 1.3) Regions file

This file contains the individual IDs present in the genealogy with the corresponding region for each invidual:

```
1	1
13	1
2	2
3	2
16	2
18	3
19	3
20	4
21	4
```
## 1.4) Regions to test file
This file contains the regions where you want to estimate the frequency of this variant. At the end of the file you have to add the region "All Probands", because this software gives you an overall estimate.

```
1
2
3
4
All Probands
```

---
#   2) Install all the dependencies required for this Software in a Python 2 environment
The software was built in Python 2 and requires some specific versions of some python packages. The original requirements for running this software can be found in the requirements.txt file
```
matplotlib==2.2.2
profilehooks==1.10.0
pandas==0.22.0
scipy==1.2.0
seaborn==0.8.1
tables==3.4.2
attrs==17.4.0
numpy==1.14
pytest==3.5.0
ipython==5.4.1
attr==0.3.1
subprocess32==3.5.3
```
However, some versions of some packages are not available right now because they are too old. In my experience, ISgen can run just fine with the following versions:

```
matplotlib==2.2.5+
profilehooks==1.10.0
pandas==0.24.2+
scipy==1.2.3+
seaborn==0.8.1+
tables==3.5.2+
attrs==17.4.0+
numpy==1.16.6+
pytest==4.6.11+
ipython==5.4.1
attr==0.3.1
subprocess32==3.5.2+
```

# 3) Running ISGen
To run all the steps of the software, you can use the script run.py, the following arguments are needed:


bash filteringagain.sh

With this script you will create a folder named "ARG_files" and the .haps and .samples files will be stored there.
