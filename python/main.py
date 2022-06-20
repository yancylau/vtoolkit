#!/usr/bin/env python

from Bio import SeqIO
import csv
from functions import nt2aa
from get_subregions import getSubregions
from determine_type import determineType




def parseSubregions(in_fasta, out_csv, amino_acid=False):
    '''
    Main function
    Arguments:
      in_fasta: FASTA file with NT/AA sequence(s)
      out_csv : CSV file to save output
    Returns:
      CSV file with parsed results
    '''
    # Fields
    fields = ['id', 'strand', 'stop_codon', 'type', 'nt', 'aa',  'v_aa', 'v_imgt_aa', 
              'fr1', 'fr2', 'cdr2', 'fr3', 'cdr1', 'cdr3', 'fr4',
              'fr1_imgt', 'cdr1_imgt', 'fr2_imgt', 'cdr2_imgt', 'fr3_imgt', 'cdr3_imgt', 'fr4_imgt', 
              'fr1_start', 'fr1_end', 'cdr1_start', 'cdr1_end', 
              'fr2_start', 'fr2_end',  'cdr2_start', 'cdr2_end', 
              'fr3_start', 'fr3_end', 'cdr3_start', 'cdr3_end', 
              'fr4_start', 'fr4_end']

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
                    subregions = getSubregions(db['aa'])
                    db.update(subregions)
            else:
                db['nt'] = record.seq
                # Detect subregions from the forward/reverse_complement sequence
                aa_fwd = nt2aa(record.seq)
                aa_rev = nt2aa(record.seq.reverse_complement())

                # TODO: Maybe only select longest one, and diff with 10 aa.
                subregions_fwd = getSubregions(aa_fwd)
                subregions_rev = getSubregions(aa_rev)
                # Determine DNA sequence strand
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
            
            


#### Usage: 
# cd ./src
# python main.py ../data/nt.fasta ../results/subregions.csv
#if __name__ == "__main__":
    #nt_fasta = "../data/nt.fasta"
    #out_csv = "../results/out.csv"
    #nt_fasta = "../gdna/collapse_unique.fasta"
    #out_csv = "../results/gdna.csv"
    #nt_fasta = sys.argv[1]
    #out_csv = sys.argv[2]
    #parseRegions(nt_fasta, out_csv)
# aa_fasta = "../data/atleast3_aa.fa"
# out_csv = "../results/atleast3_fr4.csv"


### Test on windows
# aa_fasta = "../aa_functional.fa"
# out_csv = "../subregions.csv"
# parseSubregions(aa_fasta, out_csv, amino_acid=True)


# Germline
#nt_fasta = "../data/otu.fa"
#out_csv = "../results/otu.csv"

# Rearranged
nt_fasta = "../data/cdna_atleast3.fa"
out_csv = "../results/cdna.csv"
parseSubregions(nt_fasta, out_csv, amino_acid=False)


# ### Test on linux
# aa_fasta = "../collapse_unique/aa_functional.fa"
# out_csv = "../vdj/test.csv"
# parseSubregions(aa_fasta, out_csv, amino_acid=True)