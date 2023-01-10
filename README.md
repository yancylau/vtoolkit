
# Toolkit for the annotation of antibody variable region (vtoolkit)

## Introduction

Vtoolkit determine variable (V) regions from query sequences firstly and then make numbering for V regions according to the IMGT unique numbering.  Based on the hallmarks sites located on FR2 sequence, Vtoolkit destinguish the HCAbs from conventional antibodis automaticlly.


## Functions

- Find V regions, and their subrgions, including FRs and CDRs
- IMGT numbering
- Determine V gene types



## How vtoolkits works

![figure1](figures/figure1.png)


## Usage

### Get variable region from protein sequences

Use vtoolkit to identify variable regions from query sequences.

```shell
# Usage: 
python/main.py \
  examples/aa.fa \
  examples/result.csv 

```

### Camel antibody classification: seperate heavy antibods from conventional antibodies

```shell
# Usage: 
python/main.py \
  examples/aa.fa \
  examples/result.csv 
```
