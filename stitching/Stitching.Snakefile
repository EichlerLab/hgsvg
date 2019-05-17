import os
import tempfile
import subprocess

#
# Init
#

# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)

cwd=os.getcwd()

# Tempdir
SSD_TMP_DIR = "/data/scratch/ssd"

if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
elif os.path.exists(SSD_TMP_DIR):
    TMPDIR = SSD_TMP_DIR
else:
    TMPDIR = tempfile.gettempdir()

# Config
configfile: "phasedsv.json"

with open("grid.json") as grid_file:
    gridOpts = json.load(grid_file)

for key in gridOpts:
    config[key] = gridOpts[key]

def GetFastq(fofn):
    f=open(fofn)
    return " ".join([line.strip() for line in f])

config["fastq_1_fofn"] = "NA"
config["fastq_2_fofn"] = "NA"
config["overlapsPerJob"]= 20
config["recall_bin"] = 100


# Set LD_LIBRARY_PATH

if 'ld_path' in config:
    os.environ['LD_LIBRARY_PATH'] = config['ld_path']


# Params
SLOP_FOR_SV_SEQUENCE_POSITIONS = 5000


# Chromosomes
with open(config['ref']+".fai") as faiFile:
    chroms = [l.split()[0].rstrip() for l in faiFile]


# Snakemake params
shell.prefix("source ./config.sh; ")

localrules: all


# Definitions

haps=["h0","h1"]
shortHaps=["0", "1"]

ovpIdx=range(0,config["overlapsPerJob"])


#
# All
#

rule all:
    input:
        alnFofn     = "alignments.fofn",
        alnSam      = expand("alignments.{hap}.sam", hap=haps),
        alnBam      = expand("alignments.{hap}.bam", hap=haps),
        asmBed      = expand("alignments.{hap}.bam.bed", hap=haps),
        contigBed   = expand("overlaps/overlaps.{hap}.{chrom}.ctg0.bed",hap=haps,chrom=chroms),
        asmFasta    = expand("alignments.{hap}.bam.fasta", hap=haps),
        asmOverlaps = expand("overlaps/overlap.{hap}.{chrom}.txt", hap=haps, chrom=chroms),
        asmGraphs   = expand("overlaps/overlap.{hap}.{chrom}.txt.gml", hap=haps, chrom=chroms),
        asmPaths    = expand("overlaps/overlap.{hap}.{chrom}.txt.path", hap=haps, chrom=chroms),
        asmContigs  = expand("contigs/patched.{hap}.{chrom}.fasta", hap=haps, chrom=chroms),
        alnContigs  = expand("contigs/patched.{hap}.{chrom}.fasta.sam", hap=haps, chrom=chroms),
        chrFasta    = expand("contigs.{hap}.fasta", hap=haps),
        chrFastaFai = expand("contigs.{hap}.fasta.fai", hap=haps),
        chrAln      = expand("contigs.{hap}.fasta.sam", hap=haps),
        chrBed      = expand("contigs.{hap}.fasta.sam.bed", hap=haps),
        chrBed6     = expand("contigs.{hap}.fasta.sam.bed6", hap=haps),
        chrBB       = expand("contigs.{hap}.fasta.sam.bb", hap=haps),
        gaps        = expand("stitching_hap_gaps/hap{hap}/gaps.bed", hap=shortHaps),
        indels      = expand("stitching_hap_gaps/hap{hap}/indels.orig.bed", hap=shortHaps),
        indelVCF    = expand("stitching_hap_gaps/hap{hap}/indels.orig.vcf",hap=shortHaps),
        normVCF     = expand("stitching_hap_gaps/hap{hap}/indels.norm.vcf",hap=shortHaps),
        normBed     = expand("stitching_hap_gaps/hap{hap}/indels.norm.bed",hap=shortHaps)


#
# Uncategorized rules
#

rule PrintGaps:
    input:
        contigSam="contigs.h{hap}.fasta.sam"
    output:
        gapBed="stitching_hap_gaps/hap{hap}/gaps.bed"
    params:
        sd=SD,
        ref=config["ref"],
        grid_opts=config["grid_small"]
    shell:
        """mkdir -p stitching_hap_gaps/hap{wildcards.hap}; """
        """{params.sd}/../sv/utils/PrintGaps.py {params.ref} {input.contigSam} --condense 20 --minLength 50 --outFile {output.gapBed}"""

rule ApplyPolish:
    input:
        assembly="{sample}.{hap}.fasta",
        vcf="wgsvcf_{hap}.vcf"
    output:
        polish="{sample}.{hap}.polish.fasta",
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:
        """{params.sd}/ApplyVCFPatch.py --genome {input.assembly} --vcf {input.vcf} --out {output.polish}"""

rule MakeDummyChroms:
    input:
        fasta="contigs.h0.fasta",
    output:
        dummy=expand("dummy/contig.h0.{contig}.txt",contig=chroms)
    params:
        grid_opts=config["grid_small"]
    shell:
        """samtools faidx {input.fasta}; """
        """mkdir -p dummy;\n"""
        """for ctg in `cut -f 1 {input.fasta}.fai`; do """
            """touch dummy/contig.h0.$ctg.txt; """
        """done || true"""

rule MakeDummyChromsH1:
    input:
        fasta="contigs.h1.fasta",
    output:
        dummy=expand("dummy/contig.h1.{contig}.txt",contig=chroms),
    params:
        grid_opts=config["grid_small"]
    shell:
        """samtools faidx {input.fasta}; """
        """mkdir -p dummy;\n"""
        """for ctg in `cut -f 1 {input.fasta}.fai`; do """
            """touch dummy/contig.h1.$ctg.txt; """
        """done || true"""

rule CallVariantsOnChromH0:
    input:
        fasta="contigs.h0.fasta",
        bam="index/contigs.h0.fasta.bam",
        dummy_contig="dummy/contig.h0.{contig}.txt"
    output:
        vcf="vcfs/contigs.h0.{contig}.vcf",
    params:
        grid_opts=config["grid_small"], 
    shell:
        """mkdir -p vcfs; """
        """freebayes -f {input.fasta}  -b {input.bam} --region {wildcards.contig} -m 10 --min-coverage 10 --max-coverage 40 -v {output.vcf}"""

rule CallVariantsOnChromH1:
    input:
        fasta="contigs.h1.fasta",
        bam="index/contigs.h1.fasta.bam",
        dummy_contig="dummy/contig.h1.{contig}.txt"
    output:
        vcf="vcfs/contigs.h1.{contig}.vcf",
    params:
        grid_opts=config["grid_small"], 
    shell:
        """mkdir -p vcfs; """
        """freebayes -f {input.fasta} -b {input.bam} --region {wildcards.contig} -m 10 --min-coverage 10 --max-coverage 40 -v {output.vcf}"""

rule CombineVariantsH0:
    input:
       invcf=expand("vcfs/contigs.h0.{contig}.vcf",contig=chroms),
    output:
       vcf="index/contigs.h0.fasta.vcf",
    params:
       grid_opts=config["grid_small"]
    shell:
        """grep "^#" {input.invcf[0]} > {output.vcf}; """
        """cat {input.invcf} | grep -v "^#" | sort -k1,1 -k2,2n >> {output.vcf}"""

rule CombineVariantsH1:
    input:
       invcf=expand("vcfs/contigs.h1.{contig}.vcf",contig=chroms),
    output:
       vcf="index/contigs.h1.fasta.vcf",
    params:
       grid_opts=config["grid_small"]
    shell:
        """grep "^#" {input.invcf[0]} > {output.vcf}; """
        """cat {input.invcf} | grep -v "^#" | sort -k1,1 -k2,2n >> {output.vcf}"""

rule MapHTSReadsToContigs:
    input:
       asm="contigs.h{hap}.fasta",
       bwt="index/contigs.h{hap}.fasta.bwt",
       fq1=config["fastq_1_fofn"],
       fq2=config["fastq_2_fofn"],
    output:
       bam=protected("index/contigs.h{hap}.fasta.bam"),
    params:
       grid_opts=config["grid_manycore"],
#       fq1=GetFastq(config["fastq_1_fofn"]),
#       fq2=GetFastq(config["fastq_2_fofn"])
    shell:""

#bwa mem -t 10 index/{input.asm}  '<zcat {params.fq1}' '<zcat {params.fq2}' | \
#  samtools view -bS - | \
#  samtools sort -m2G -@ 4 -T $TMPDIR/aln -o {output.bam}


rule MakeMemIndex:
    input:
       asm="contigs.h{hap}.fasta"
    output:
       bwt="index/contigs.h{hap}.fasta.bwt"
    params:
       grid_opts=config["grid_large"]
    shell:
        """mkdir -p index; """
        """cd index && ln -s ../{input.asm} . && cd .. ;"""
        """bwa index index/{input.asm}"""

rule MakeIndels:
    input:
       asmSam="contigs.h{hap}.fasta.sam",
    output:
       indels="stitching_hap_gaps/hap{hap}/indels.orig.bed"
    params:
       grid_opts=config["grid_small"],
       sample=config["sample"],
       sd=SD,
       ref=config["ref"],
    shell:
        """mkdir -p stitching_hap_gaps/hap{wildcards.hap}; """
        """{params.sd}/../sv/utils/PrintGaps.py {params.ref} {input.asmSam} --minLength 2 --maxLength 49 --ignoreHP 5 --outFile {output.indels}"""

rule ConvertIndelBedToVCF:
    input:
        indelBed="stitching_hap_gaps/hap{hap}/indels.orig.bed"
    output:
        indelVCF="stitching_hap_gaps/hap{hap}/indels.orig.vcf"
    params:
        grid_opts=config["grid_long"],
        ref=config["ref"],
        sample=config["sample"]
    shell:
        """{SD}/../sv/utils/variants_bed_to_vcf.py --bed {input.indelBed} --ref {params.ref} --sample {params.sample} --type indel --vcf /dev/stdout | """
        """bedtools sort -header """
        """> {output.indelVCF}"""

rule NormIndelVCF:
    input:
        indelVCF="stitching_hap_gaps/hap{hap}/indels.orig.vcf"
    output:
        indelNormVCF="stitching_hap_gaps/hap{hap}/indels.norm.vcf"
    params:
        grid_opts=config["grid_small"],
        ref=config["ref"],
    shell:
        """vt normalize -n -r {params.ref} -o {output.indelNormVCF} {input.indelVCF}"""

rule NormIndelVCFToBed:
    input:
        indelNormVCF="stitching_hap_gaps/hap{hap}/indels.norm.vcf"
    output:
        indelNormBed="stitching_hap_gaps/hap{hap}/indels.norm.bed"
    params:
        grid_opts=config["grid_small"],
        ref=config["ref"],
    shell:
        """{SD}/../sv/utils/variants_vcf_to_bed.py --vcf {input.indelNormVCF} --out {output.indelNormBed}"""

rule CombineHapIndels:
    input:
        indelNormBed=expand("stitching_hap_gaps/hap{hap}/indels.norm.bed",hap=shortHaps)
    output:
        indelBed="stitching_hap_gaps/diploid/indels.bed"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:
        """{params.sd}/../sv/utils/MergeHaplotypes.sh {input.indelNormBed} {output.indelBed} "svType svLen svSeq" 0.8"""
    
rule MakeAnnotation:
    input:
        asmSam=expand("contigs.{hap}.fasta.sam",hap=haps)
    output:
        annotation="stitching_hap_gaps/diploid/insertions.bed"
    params:
        grid_opts=config["grid_small"]
    shell:
        "make -f " + SD + "/DiploidAnnotation.mak H0SAM=contigs.h0.fasta.sam H1SAM=contigs.h1.fasta.sam DIR=stitching_hap_gaps -j 2"
    
rule MakeChrAsmFasta:
    input:
        asmContig=expand("contigs/patched.{{hap}}.{chrom}.fasta", chrom=chroms)
    output:
        asmFasta="contigs.{hap}.fasta"
    params:
        grid_opts=config["grid_small"],
    shell:
        "cat {input.asmContig} > {output.asmFasta}"

rule MakeChrAsmFastaFai:
    input:
        asmFasta="contigs.{hap}.fasta"
    output:
        asmFastaFai="contigs.{hap}.fasta.fai"
    params:
        grid_opts=config["grid_small"],
    shell:
        "samtools faidx {input.asmFasta}"
    
rule MakeAsmAln:
    input:
        asmContigs=expand("contigs/patched.{{hap}}.{chrom}.fasta.sam", chrom=chroms),
        aln="alignments.{hap}.bam"
    output:
        asmSam="contigs.{hap}.fasta.sam"
    params:
        grid_opts=config["grid_small"],
    run:

        with open(output.asmSam, 'w') as out_file:

            # Fetch and write header
            samtools_stdout = subprocess.run(
                [
                    'samtools',
                    'view',
                    '-H',
                    'alignments.h0.bam'
                ],
                stdout=subprocess.PIPE
            ).stdout.decode()

            for line in samtools_stdout.split('\n'):
                line = line.strip()

                if not line:
                    continue

                out_file.write(line)
                out_file.write('\n')

            # Iterate records and write
            for asm_sam_file_name in input.asmContigs:
                with open(asm_sam_file_name, 'r') as in_file:
                    for line in in_file:
                        line = line.strip()

                        if not line:
                            continue

                        if line.startswith('@'):
                            continue

                        out_file.write(line)
                        out_file.write('\n')


rule MakeContigAsmAln:
    input:
        asmFasta="contigs/patched.{hap}.{chrom}.fasta"
    output:
        asmSam="contigs/patched.{hap}.{chrom}.fasta.sam"
    params:
        grid_opts=config["grid_quad"],
        ref=config['ref'],
        sd=SD,
        td=TMPDIR        
    shell:
        """{params.sd}/MapContigs.py --contigs {input.asmFasta} --ref {params.ref} --tmpdir $TMPDIR --blasr {params.sd}/../../dep/bin/blasrmc --out {output.asmSam} --nproc 4"""

rule MakeChrAsmBed:
    input:
        asmSam="contigs.{hap}.fasta.sam"
    output:
        asmBed="contigs.{hap}.fasta.sam.bed"
    params:
        grid_opts=config["grid_small"],
        hgsvg=SD+ "/.."
    shell:
        "{params.hgsvg}/mcutils/src/samToBed {input.asmSam}  --reportIdentity | bedtools sort  > {output.asmBed}"

rule MakeAsmBB:
    input:
        asmBed="contigs.{hap}.fasta.sam.bed6"
    output:
        asmBB="contigs.{hap}.fasta.sam.bb"
    params:
        grid_opts=config["grid_small"],
        ref=config['ref']
    shell:
        "bedToBigBed {input.asmBed} {params.ref}.fai {output.asmBB} -type=bed6"

rule MakeChrAsmBed6:
    input:
        asmBed="contigs.{hap}.fasta.sam.bed"
    output:
        asmBed6="contigs.{hap}.fasta.sam.bed6"
    params:
        grid_opts=config["grid_small"],
        hgsvg=SD+ "/.."
    shell:
        "{params.hgsvg}/utils/tracks/SamBedToBed6.py {input.asmBed} {output.asmBed6} "


#
# Overlaps
#

rule MakeAsmContigs:
    input:
        asmPath="overlaps/overlap.{hap}.{chrom}.txt.path",
        asmOverlap="overlaps/overlap.{hap}.{chrom}.txt",
        asmFasta="alignments.{hap}.bam.fasta"
    output:
        asmContig="contigs/patched.{hap}.{chrom}.fasta"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:
        """{params.sd}/PatchPaths.py {input.asmOverlap} {input.asmFasta} {input.asmPath} {output.asmContig}"""

rule MakeAsmPaths:
    input:
        asmOverlap="overlaps/overlap.{hap}.{chrom}.txt",
        asmOverlapGraph="overlaps/overlap.{hap}.{chrom}.txt.gml"
    output:
        asmPath="overlaps/overlap.{hap}.{chrom}.txt.path"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:
        """{params.sd}/OverlapGraphToPaths.py {input.asmOverlap} {input.asmOverlapGraph} {output.asmPath}"""

rule MakeAsmGraphs:
    input:
        asmOverlap="overlaps/overlap.{hap}.{chrom}.txt"
    output:
        asmOverlapGraph="overlaps/overlap.{hap}.{chrom}.txt.gml"
    params:
        grid_opts=config["grid_small"],
        sd=SD
    shell:
        """{params.sd}/OverlapsToGraph.py {input.asmOverlap} --out {output.asmOverlapGraph}"""

#subworkflow AsmOverlapsWorkflow:
#    snakefile: SD +"/MakeAsmOverlaps.Snakefile"
#    workdir:cwd
#

rule MakeAsmOverlaps:
    input:
        bed="overlaps/overlaps.{hap}.{chrom}.ctg0.bed",
        asm="alignments.{hap}.bam.fasta",
        splitAsmOverlaps=expand("overlaps/split_{{chrom}}/overlaps.{{hap}}.{{chrom}}.ctg0.{start}.txt",start=ovpIdx)
    output:
        asmOverlap="overlaps/overlap.{hap}.{chrom}.txt"
    params:
        grid_opts=config["grid_small"]
    shell:
        "cat {input.splitAsmOverlaps} > {output.asmOverlap}"

rule MakeSplitOverlaps:
    input:
        bed="overlaps/split_{chrom}/overlaps.{hap}.{chrom}.ctg0.{start}.bed",
        asm="alignments.{hap}.bam.fasta"
    output:
        splitAsmOverlaps="overlaps/split_{chrom}/overlaps.{hap}.{chrom}.ctg0.{start}.txt"
    params:
        grid_opts=config["grid_manycore"],
        ovps=config["overlapsPerJob"],
        sd=SD
    shell:
        """mkdir -p $TMPDIR; """
        """mkdir -p overlaps/split_{wildcards.chrom}; """
        """{params.sd}/OverlapContigsOrderedByBed.py {input.bed} {input.asm} """
            """--chrom {wildcards.chrom} """
            """--out {output.splitAsmOverlaps} """
            """--nproc 12 """
            """--tmpdir $TMPDIR """
            """--blasr {params.sd}/../../dep/bin/blasrmc """
            """--path {params.sd}"""

rule SplitAsmOverlaps:
    input:
        bed="overlaps/overlaps.{hap}.{chrom}.ctg0.bed",
    output:
        split=expand("overlaps/split_{{chrom}}/overlaps.{{hap}}.{{chrom}}.ctg0.{start}.bed", start=ovpIdx),
    params:
        grid_opts=config["grid_small"],
        sd=SD,
        nOvp=config['overlapsPerJob']
    shell:
        """mkdir -p overlaps/split_{wildcards.chrom}; """
        """{params.sd}/SplitBedFile.py --bed {input.bed} --n {params.nOvp} --overlap 5 --base overlaps/split_{wildcards.chrom}/overlaps.{wildcards.hap}.{wildcards.chrom}.ctg0"""

rule MakeContigBed:
    input:
        asmBed = "alignments.{hap}.bam.bed"
    output:
        contigBed = expand("overlaps/overlaps.{{hap}}.{chrom}.ctg0.bed",chrom=chroms)
    params:
        grid_opts=config["grid_small"],
        ch=chroms
    shell:
        """for c in {params.ch} ; do \n"""
            """    egrep "^$c\t" {input.asmBed} | """
            """ egrep "/0" > overlaps/overlaps.{wildcards.hap}.$c.ctg0.bed;\n"""
        """done || true"""

rule MakeAsmBed:
    input:
        asmBam = "alignments.{hap}.bam"
    output:
        asmBed = "alignments.{hap}.bam.bed"
    params:
        grid_opts=config["grid_small"]
    shell:
        """samtools view {input.asmBam} | """
        """{SD}/../../dep/bin/samToBed /dev/stdin --reportIdentity """
        """> {output.asmBed}"""

rule MakeAsmFasta:
    input:
        asmBam = "alignments.{hap}.bam"
    output:
        asmFasta = "alignments.{hap}.bam.fasta"
    params:
        grid_opts=config["grid_small"]
    shell:
        """samtools view {input.asmBam} | """
        """awk '{{ print \">\"$1; print $10;}}' | """
        """fold | """
        """sed '/^$/d' """
        """> {output.asmFasta}; """
        """samtools faidx {output.asmFasta}"""


#
# Init alignments
#

rule MakeAlnBam:
    input:
        sam="alignments.{hap}.sam",
    output:
        bam="alignments.{hap}.bam",
    params:
        sd=SD,
        grid_opts=config["grid_small"]
    shell:
        """samtools view -bS {input.sam} | samtools sort -T $TMPDIR/tmp.{wildcards.hap}.sam -o {output.bam}; """
        """samtools index {output.bam}"""

rule MakeAlnSam:
    input:
        fofn="alignments.fofn"
    output:
        sams=expand("alignments.{hap}.sam",hap=haps)
    params:
        sd=SD,
        grid_opts=config["grid_small"]
    shell:
        """{params.sd}/FilterSamByHaplotype.py {input.fofn} --h0 {output.sams[0]} --h1 {output.sams[1]} --header {params.sd}/header.sam"""

rule MakeFofn:
    output:
        fofn="alignments.fofn"
    params:
        grid_opts=config["grid_small"]
    shell:
        """ls samfiles | awk '{{ print "samfiles/"$1;}}' > {output.fofn}"""
