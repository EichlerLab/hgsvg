#!/usr/bin/env python
import argparse
import datetime
import numpy as np
import pandas as pd
import pysam
import sys
import math

def calculate_variant_quality(variant):
    try:
        return int(min(100, round(-10 * np.log10(1 - (variant.contig_support / float(variant.contig_depth))) * np.log(variant.contig_depth), 0)))
    except ZeroDivisionError:
        return 0

def GetGenotype(gt):
    if gt == "HAP0":
        return "1|0"
    elif gt == "HAP1":
        return "0|1"
    elif gt == "HOM":
        return "1|1"
    elif gt == "HET":
        # unphased
        return "1/0"
    else:
        # assume unphased homozygous
        return "1/1"


def GetStart(row):
    if row["svType"] == "insertion":
        return row["tStart"]
    else:
        # deletion, gives last ref base before event.
        return row["tStart"]

def GetType(val):
    if len(val)> 3:
        return val[0:3].upper()
    else:
        return "NONE"

def GetPass(val):
    if val > 3:
        return "PASS"
    else:
        return "FAIL"

def GetSeq(val):
    if type(val) == float:
        return "NA"
    else:
        return val

def GetAltSeq(row, genome):
    seq = GetSeq(row["svSeq"]).upper()
    if row["svType"] == "insertion":
        refPrefix = genome.fetch(row["#chrom"], row["origTStart"]-1, row["origTStart"]).upper()
        return refPrefix + seq.upper()
    else:
        refPrefix = genome.fetch(row["#chrom"], row["origTStart"]-1, row["origTStart"]).upper()        
        return refPrefix

def ParseSVLen(val):
    if str(val).find(",") != -1:
        return [int(i) for i in val.split(",")]
    else:
        try:
            i=abs(int(val))
        except:
            import pdb
            pdb.set_trace()
            print(val)
            
        return abs(int(val))
    
def GetSVLen(row):
    retval= None

    if row["svType"] == "locus":
        retval = row["svLen"]
    elif row["svType"] == "deletion":
        retval = -1*int(row["svLen"])
    else:
        retval = row["svLen"]

    if retval is None or retval == "":
        import pdb
        pdb.set_trace()
        print(row["svType"])
    return retval

def GetRefSeq(row, genome, oneBase=False):
    refLen = 1;
    refPos = row["origTStart"]
    if row["svType"] == "deletion":
        if oneBase == False:
            refLen = row["svLen"] + 1
    elif row["svType"] == "locus":
        refLen =  row["tEnd"] - refPos
        # deletion sequence starts at base before del event
    refPos-=1

    seq = genome.fetch(row["#chrom"], refPos, refPos + refLen).upper()
    return seq.upper()



def convert_bed_to_vcf(bed_filename, reference_filename, vcf_filename, sample, variant_type):
    # Get variants.
    if variant_type == "sv":
        columns = (0, 1, 2, 3, 4, 5, 7, 9, 10, 12, 14)
        names = ("chrom", "start", "end", "sv_call", "event_size", "sv_sequence", "contig", "contig_start", "contig_end", "genotype", "repeat_type")
        fmt     = ["GT"]
    elif variant_type == "indel":
        # chr1    94824   94827   3       Cttttcttttttttt 1       1       29.04   deletion
        columns = (0, 1, 2, 3, 4, 5, 6, 7, 8)
        names = ("chrom", "start", "end", "event_size", "sv_sequence", "contig_support", "contig_depth", "depth", "sv_call")
        fmt=["GT"]
    elif variant_type == "inversion":
        columns = (0, 1, 2, 3, 4, 5)
        names = ("chrom", "start", "end", "sv_call", "contig_support", "contig_depth")
        fmt=[""]
    else:
        raise Exception("Unsupported variant type: %s" % variant_type)


    
    calls = pd.read_table(bed_filename,low_memory=False, keep_default_na=False, index_col=False, header=0)#, header=None, usecols=columns, names=names)
    
    calls["sample_name"] = sample
    calls["call_id"] = "."
    calls["quality"] = "30" #calls.apply(calculate_variant_quality, axis=1)
    calls["filter"] = "PASS"


    # Make sure the sv length and sv sequence agree

#    calls["svLen"] = calls.apply(lambda row: len(GetSeq(row["svSeq"])), axis=1)

    pd.to_numeric(calls["tStart"])
    pd.to_numeric(calls["tEnd"])
    # Get the reference base at the position of the variant start.
    reference = pysam.FastaFile(reference_filename)
    calls["reference"] = calls.apply(lambda row: reference.fetch(row["#chrom"], row["tStart"], row["tStart"] + 1).upper(), axis=1)
    
    # Update start position to be 1-based.
    calls["origTStart"] = calls["tStart"]
    calls["CHROM"]  = calls["#chrom"]

    calls["POS"] = calls.apply(lambda row: GetStart(row), axis=1)
    if args.addci is not None:
        calls["CIPOS"] = ["-{},{}".format(args.addci, args.addci)]*len(calls)
        calls["CIEND"] = ["-{},{}".format(args.addci, args.addci)]*len(calls)        
    # Build an INFO field for each call.
    calls["svShort"] = calls.apply(lambda row: GetType(row["svType"]), axis=1)
    

    
    if variant_type == "sv":
        infoKeys= [("END", "tEnd"),("SVTYPE", "svShort"),("SVLEN", "svLen"),("CONTIG", "qName"),("CONTIG_START", "qStart"),("CONTIG_END", "qEnd"),("SEQ", "svSeq")]
        if "is_trf" in calls:
            infoKeys.append(("IS_TRF", "is_trf"))

        if args.addci is not None:
            infoKeys.append(("CIEND", "CIEND"))
            infoKeys.append(("CIPOS","CIPOS"))
        
        if len(args.fields) > 0:
            extraKeys = [(args.fields[i], args.fields[i+1]) for i in range(0,len(args.fields),2)]
            infoKeys += extraKeys
        if args.seq:
            calls["reference"] = calls.apply(lambda row: GetRefSeq(row, reference), axis=1)
            calls["alt"] = calls.apply(lambda row: GetAltSeq(row, reference),axis=1)
        else:
            calls["reference"] = calls.apply(lambda row: GetRefSeq(row, reference, 1), axis=1)
            calls["alt"] = calls.apply(lambda row: "<%s>" % row.svType[:3].upper(), axis=1) 

        calls["svLen"] =  calls.apply(lambda row: GetSVLen(row), axis=1)
            
        calls["info"] = calls.apply(
            lambda row: ";".join(
                ["=".join(map(str, (item[0], row[item[1]])))
                 for item in ( infoKeys )]
            ),
            axis=1
        )
#        import pdb
#        pdb.set_trace()
        calls["svLen"] = calls.apply(lambda row: ParseSVLen(row["svLen"]), axis=1)
        calls["format"] = ":".join(fmt)
        if "hap" in calls:
            calls["genotype"] = calls.apply(lambda row: GetGenotype(row.hap), axis=1)
        else:
            calls["genotype"] = ["./."]* len(calls["tEnd"])
        
    elif variant_type == "indel":
        calls["reference"] = calls.apply(lambda row: GetRefSeq(row, reference), axis=1)
        calls["alt"] = calls.apply(lambda row: GetAltSeq(row, reference),axis=1)
        calls["format"] = ":".join(fmt)
        
        if "hap" in calls:
            calls["genotype"] = calls.apply(lambda row: GetGenotype(row.hap), axis=1)
        else:
            calls["genotype"] = ["./."]* len(calls["tEnd"])
        calls["svLen"] =  calls.apply(lambda row: GetSVLen(row), axis=1)        
        calls["info"] = calls.apply(
            lambda row: ";".join(
                ["=".join(map(str, item))
                 for item in (
                        ("END", row["tEnd"]),
                        ("SVTYPE", row["svType"]),
                        ("SVLEN", row["svLen"]),
                        ("SAMPLES", row["sample_name"]),
                        ("SEQ", row["svSeq"])
                 )]
            ),
            axis=1
        )
    elif variant_type == "inversion":
        calls["alt"] = "<INV>"
        calls["info"] = calls.apply(
            lambda row: ";".join(
                ["=".join(map(str, item))
                 for item in (
                        ("END", row["tEnd"]),
                        ("SVTYPE", row["svType"]),
                         ("SVLEN", row["svLen"]),
                        ("SAMPLES", row["sample_name"]),
                 )]
            ),
            axis=1
        )

    simple_calls = calls[["#chrom", "POS", "call_id", "reference", "alt", "quality", "filter", "info", "format", "genotype"]].rename({"#chrom": "#CHROM", "reference": "REF", "call_id": "ID", "quality": "QUAL", "info": "INFO", "alt": "ALT", "filter": "FILTER", "format":"FORMAT", "genotype": sample}, axis=1)

    faiFile = open(args.reference + ".fai")
    fai = []
    for line in faiFile:
        vals = line.split()
        fai.append([vals[0], vals[1]])
        
    # Save genotypes as tab-delimited file.
    with open(vcf_filename, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##fileDate=%s\n" % datetime.date.strftime(datetime.date.today(), "%Y%m%d"))
        vcf.write("##source={}\n".format(args.source))
        vcf.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' + "\n")
        vcf.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">' + "\n")
        vcf.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">' + "\n")
        vcf.write('##INFO=<ID=CONTIG,Number=1,Type=String,Description="Name of alternate assembly contig">' + "\n")
        vcf.write('##INFO=<ID=CONTIG_START,Number=1,Type=Integer,Description="Start coordinate of this variant in the alternate assembly contig">' + "\n")
        vcf.write('##INFO=<ID=CONTIG_END,Number=1,Type=Integer,Description="End coordinate of this variant in the alternate assembly contig">' + "\n")
        vcf.write('##INFO=<ID=SEQ,Number=1,Type=String,Description="Sequence associated with variant">' + "\n")
        vcf.write('##INFO=<ID=SAMPLES,Number=1,Type=String,Description="Sample associated with variant">' + "\n")
        for i in range(0,len(fai)):
            vcf.write("##contig=<ID={},length={}>".format(fai[i][0], fai[i][1]) + "\n")
        vcf.write("##SAMPLE=<ID={}>\n".format(args.sample))
        vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + "\n")
        if args.info is not None:
            vcf.write("\n".join(args.info)+"\n")
        vcf.flush()
        simple_calls.to_csv(vcf, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", help="input BED file of variant calls")
    parser.add_argument("--reference", help="FASTA file for reference used with SMRT SV", required=True)
    parser.add_argument("--vcf", help="output VCF file of variant calls")
    parser.add_argument("--seq", help="output ref and alt sequences, not <DEL> and <INS> for svs.", action='store_true', default=False)
    parser.add_argument("--sample", help="name of sample with variants")
    parser.add_argument("--type", help="variant call type", choices=("sv", "indel", "inversion"))
    parser.add_argument("--fields", help="Additional  fields to add", nargs="+", default=[])
    parser.add_argument("--info", help="Additional info descriptions to match fields.", default=None, nargs="+")    
    parser.add_argument("--source", help="Source of variants file.", default="SMRT_SV"),
    parser.add_argument("--addci", help="Add CIPOS=-v,v;CIEND=-v,v where v is the parameter value", default=None,type=int)
    args = parser.parse_args()

    convert_bed_to_vcf(args.bed, args.reference, args.vcf, args.sample, args.type)
