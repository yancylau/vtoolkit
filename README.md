# vtoolkit: Toolkit for the annotation of antibody variable (V) region
Update: January 10, 2023 (Version 1.0)


## Introduction
Vtoolkit identify V regions of antibodies firstly, and number V regions according to the IMGT unique numbering scheme. Based on the hallmark sites (42, 49, 50 and 52) on FR2, vtoolkit is able to separate variable genes of HCAbs (VHHs) from variable genes of conventional antibodies (VHs).


## Functions
- Find V regions, and their subrgions, including FRs and CDRs
- IMGT numbering
- Determine V gene types


## How vtoolkit works
<p align="center">
  <img width="600"  src="figures/figure1.png">
</p>

### V (region) detector
V detector was packed in vtoolkit to identify V regions of antibody heavy chain. It consists of 4 FR detectors, FR1/2/3/4 detector, which were fitted with Naive Bayes Classifier (scikit-learn 1.2.0). Initially, FRs were detected with FR1/2/3/4 detectors. After the detection of FRs, CDRs and the whole V region could be easily inferred from the query sequences (Figure 1). Then, vtoolkit will number V regions according to IMGT scheme.

### V gene classifier
In vtoolkit, V classifier was applied to determine the types of the V gene, which were fitted with Naive Bayes Classifier (scikit-learn 1.2.0). The hallmark sites of 42, 49, 50 and 52 were inclued as features in V calssifer.


## Usage
### Identify V regions from antibodies (protein sequences)
Use vtoolkit to identify variable regions from query sequences.
```shell
python/main.py \
  examples/refs.faa \
  examples/refs.annotated.csv 
```

### Camel V gene classification: separate VHHs from VHs
```shell
python/main.py \
  examples/camel.faa \
  examples/camel.annotated.csv 
```

## References


## Release Notes
v1.0, 2023.1.10, original version
