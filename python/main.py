#!/usr/bin/env python3

import sys
from Bio import SeqIO
import csv
from functions import nt2aa
from get_v import getV
from determine_v_type import determineType



def parseSeqs(in_fasta, out_csv, amino_acid=True):
    '''
    Main function
    Arguments:
      in_fasta: fasta file with AA/NT sequences
      out_csv : csv file to save output
    Returns:
      csv file with parsed results
    '''
    # Fields
    fields = ['id', 
              'strand', 
              'stop_codon', 
              'type', 
              'nt', 
              'aa',  
              'v_aa', 
              'v_imgt_aa', 
              'fr1', 
              'fr2', 
              'fr3', 
              'fr4',
              'cdr1', 
              'cdr2', 
              'cdr3', 
              'cdr3_length',
              'fr1_imgt', 
              'fr2_imgt', 
              'fr3_imgt', 
              'fr4_imgt', 
              'cdr1_imgt', 
              'cdr2_imgt', 
              'cdr3_imgt', 
              'fr1_start', 'fr1_end', 
              'fr2_start', 'fr2_end', 
              'fr3_start', 'fr3_end', 
              'fr4_start', 'fr4_end',
              'cdr1_start', 'cdr1_end', 
              'cdr2_start', 'cdr2_end', 
              'cdr3_start', 'cdr3_end']

    with open(out_csv, 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=',')
        writer.writeheader()

        for record in SeqIO.parse(in_fasta, "fasta"):
            db = {}
            db['id'] = record.id

            # Get subregions
            if amino_acid:
                db['aa'] = record.seq
                # Check stop codon
                if db['aa'].find("*") != -1:
                    db['stop_codon'] = "T"
                else:
                    db['stop_codon'] = "F"
                    subregions = getV(db['aa'])
                    db.update(subregions)
            else:
                db['nt'] = record.seq
                # Get variabe region from the forward/reverse_complement sequence
                aa_fwd = nt2aa(record.seq)
                aa_rev = nt2aa(record.seq.reverse_complement())

                # TODO: Maybe only select longest one, and diff with 10 aa.
                subregions_fwd = getV(aa_fwd)
                subregions_rev = getV(aa_rev)
                # Determine DNA strand
                if len(set(subregions_fwd.values())) == len(set(subregions_rev.values())):
                    # When FR not found
                    db['strand'] = ""
                    db['aa'] = ""
                    subregions = subregions_fwd
                if len(set(subregions_fwd.values())) > len(set(subregions_rev.values())): 
                    db['strand'] = "+"
                    db['aa'] = aa_fwd
                    subregions = subregions_fwd
                if len(set(subregions_fwd.values())) < len(set(subregions_rev.values())): 
                    db['strand'] = "-"
                    db['aa'] = aa_rev
                    subregions = subregions_rev
                db.update(subregions)

            # Determine sequence type
            if len(db['fr2']) == 17:
                try: 
                    db['type'] = determineType(str(db['fr2_imgt']))
                except:
                    db['type'] = "NA"

            writer.writerow(db)
            
            


## Usage: 
# cd vtoolkit
# python python/main.py test/v3_aa.fa test/v3_aa_annotation.csv

if __name__ == "__main__":
    nt_fasta = sys.argv[1]
    out_csv = sys.argv[2]
    parseSeqs(nt_fasta, out_csv)

## V3 genes
# aa_fasta = "test/v3_aa.fa"
# out_csv = "test/v3_aa_annotation.csv"
# parseSeqs(aa_fasta, out_csv, amino_acid=True)

## Camel genes
# aa_fasta = "test/camel_aa.fa"
# out_csv = "test/camel_aa_annotation.csv"
# parseSeqs(aa_fasta, out_csv, amino_acid=True)




### Test on windows
# aa_fasta = "../aa_functional.fa"
# out_csv = "../subregions.csv"
# parseSubregions(aa_fasta, out_csv, amino_acid=True)

# ### Test on linux
# aa_fasta = "../collapse_unique/aa_functional.fa"
# out_csv = "../vdj/test.csv"
# parseSubregions(aa_fasta, out_csv, amino_acid=True)
