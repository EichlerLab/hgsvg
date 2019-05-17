import os
import tempfile

#
# A little complicated to find the temp dir
#
SSD_TMP_DIR = "/data/scratch/ssd"
if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
elif os.path.exists(SSD_TMP_DIR):
    TMPDIR = SSD_TMP_DIR
else:
    TMPDIR = tempfile.gettempdir()

configfile: "merge.json"


SLOP_FOR_SV_SEQUENCE_POSITIONS = 5000

faiFile = open(config['ref']+".fai")
chroms = [l.split()[0].rstrip() for l in faiFile]

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
SD=SNAKEMAKE_DIR
cwd=os.getcwd()

shell.prefix(". {SNAKEMAKE_DIR}/config.sh; ")

haps=["h0","h1"]
parents=["fa", "mo"]
mssmHaps=["1", "2"]

allMSSM = [config["mssm-h0"],config["mssm-h1"],config["mssm-dn"]]
mssmhap = [config["mssm-h0"],config["mssm-h1"]]

fullMSSM={ os.path.basename(f) : f for f in allMSSM}
allMSSMBase = {f: os.path.basename(f) for f in allMSSM}
mssmBaseToPath = {os.path.basename(f): f for f in allMSSM}
sources=["uw","mssm"]
longOps=["deletion","insertion"]
shortOps=["del","ins"]
filtered=["uw.bed.filt", "mssm.bed.filt"]
#shortOpToBnOp = { "del": "DEL", "ins" : "INS" }
shortToLong= { "del" : "deletion", "ins" : "insertion" }

gapdir=config["gapdir"]
rule all:
    input:
        mssmLocal  = expand("local.{mssm}", mssm=allMSSMBase.values()),
        mssmUW     = expand("{mssm}.uw",mssm=allMSSMBase.values()),
        mssmSort   = expand("{mssm}.sorted",mssm=allMSSMBase.values()),
        mssmSortDD = expand("{mssm}.sorted.dedup",mssm=allMSSMBase.values()),
        mssmSortBN = expand("{mssm}.sorted.bntab",mssm=allMSSMBase.values()),        
        mssmSortClusters   = expand("{mssm}.sorted.window_clusters",mssm=allMSSMBase.values()),        
        mssmBnSupport=expand("bn_support.{mssm}.tab",mssm=allMSSMBase.values()),
        mssmBPSupport=expand("{mssm}.uw.cov",mssm=allMSSMBase.values()),
        mssmToUWFilt = expand("{mssm}.uw.bed.filt",mssm=allMSSMBase.values()),
        mssmBNAnnotation   = expand("{mssm}.uw.bed.bn",mssm=allMSSMBase.values()),
        bnCalls="calls_bn.bed",
        bnCallsOp=expand("calls_bn.{op}.bed",op=longOps),
        svBnOvp=expand("bn_overlap.{caller}.bed.filt.{op}.bed",caller=sources, op=shortOps),
        svFilt=expand("bn_filt.{source}.{op}.bed",source=sources, op=shortOps),
        optSvCalls=expand("opt_bn_overlap.bed.filt.{op}.bed",op=shortOps),
        uwBNAnnotation = "uw.bed.bn",
        mergedMSSM = "mssm.bed.filt",        
        uwFilt="uw.bed.filt",
        simpleReciprocal=expand("reciprocal.{op}.bed",op=shortOps),
        simpleWindow=expand("window.{op}.bed",op=shortOps),        
        svbb=expand("{source}.{op}.bb", source=sources, op=longOps),
        svs=expand("{filt}.{op}.bed", filt=filtered,op=shortOps),
        svsbn=expand("{filt}.{op}.bed.bn", filt=filtered,op=shortOps),
        overlaps=expand("{filt}.{op}.overlaps.bed", filt=filtered,op=shortOps),
        svAnnot=["uw.bed.filt.ins.mssm.bed.filt-bn.bed","uw.bed.filt.del.mssm.bed.filt-bn.bed", "mssm.bed.filt.ins.uw.bed.filt-bn.bed","mssm.bed.filt.del.uw.bed.filt-bn.bed"],
        svDist=expand("mssm.bed.filt.{op}.bed.uw-dist",op=shortOps),
        svUWDist=expand("uw.bed.filt.{op}.bed.mssm-dist",op=shortOps),
        mergedCalls=expand("merged.{op}.bed",op=shortOps),
        mergedPreInvBed="sv_calls.pre-inv.bed",
        clusterCount=[config["uwsv"]+".cluster-count"] + expand("{mssm}.uw.cluster-count",mssm=allMSSMBase.values()),
        clusterSummary="cluster-summary.txt",
        annotatedMergedBed="sv_calls.bed.annotated",
        annotatedMergedBedFixed="sv_calls.bed.annotated.tr_fixed",        
        mergedBedWithLoci="sv_calls.bed.annotated.tr_fixed.loci",
#        excludeSimpleAndMEI="excluding_simple_and_mei.bed",
#        UniqueExons="excluding_simple_and_mei_exon.bed",        
        svMSSMSource="sv_calls.base.bed",
        clusteredSVs="cluster-count.tsv",
        svUWSource="sv_calls.uw-source.bed",
        svPreSourceCalls="sv_calls.pre-source.bed",
        mergedVCF="sv_calls.vcf",
        mergedVCFgz=expand("{sample}.sv_calls.vcf.gz", sample=config["sample"]),
        mergedBB=expand("{sample}.{op}.bb",sample=config["sample"],op=shortOps),
        inv=expand("{dir}/inversions/inversions.{hap}.bed", dir=config["alnDir"], hap=haps),
        invMerged=expand("inversions.merged.{hap}.bed", hap=haps),
        invHaps=expand("{sample}.inversions.UW.bed",sample=config["sample"]),
        paCov=expand("parent.cov.{pa}.bed", pa=parents),
        isectCounts=expand("{op}.intersect.txt",op=shortOps),
        merged=expand("merged.sources.{op}.txt",op=shortOps),
        bnCompare=["sv_calls.bed.ins.bn.bed", "sv_calls.bed.ins.bn.bed", "sv_calls.bed.by-pb.bn.bed"],
        bnPlot=expand("BioNano.PacBio.Agreement.{sample}.pdf",sample=config["sample"]),
        plot=config["sample"]+".parental_coverage.pdf",
        plotnotrf=config["sample"]+".parental_coverage.no_trf.pdf",        
        hethom="het_hom_summary.txt",
        alt="sv_calls.bed.annotated.alt",
        reportMSSM=expand("{mssm}.uw.bed.filt.report", mssm=allMSSMBase.values()),
        optSummary=expand("opt_bn_overlap.bed.filt.{op}.bed.summary",op=shortOps),
        exons=expand("{sample}.exons.{op}.bed",sample=config["sample"],op=longOps),
#        sv_intv="sv_calls.intv",
#        sv_clust="sv_calls.clust",
#        sv_tr="sv_calls.clust.tr",
#        sv_tr_annot="sv_calls.clust.tr.annot",
#        sv_tr_bg="sv_calls.clust.tr.bg",        
#        wide="sv_calls.clust.wide",
        comb=expand("{mssm}.hap{hap}.bed.uw.bed.filt.comb",mssm=config["sample"],hap=mssmHaps),        
        trNet=expand(gapdir+"/hap{hap}/tr_net.tab",hap=mssmHaps),
        trZyg=expand(gapdir+"/hap{hap}/tr_net.tab.zyg",hap=mssmHaps),
        trZygBed=expand(gapdir+"/hap{hap}/tr_net.tab.zyg.bed",hap=mssmHaps),
        hapInDir=expand(gapdir+"/hap{hap}/gaps.recalled",hap=mssmHaps),
        hapInDirFilt=expand(gapdir+"/hap{hap}/gaps.recalled.filt",hap=mssmHaps),        
        filtTR=expand(gapdir+"/hap{hap}/gaps.recalled.noclust.tr_bed",hap=mssmHaps),
        notr=expand(gapdir+"/hap{hap}/gaps.recalled.noclust.notr",hap=mssmHaps),
        hettr=expand(gapdir+"/hap{hap}/gaps.recalled.noclust.hettr",hap=mssmHaps),
        gapsClust=expand(gapdir+"/hap{hap}/gaps.recalled.clust",hap=mssmHaps),
        gapsNoClust=expand(gapdir+"/hap{hap}/gaps.recalled.noclust",hap=mssmHaps),
        meiReport="mei_report.tab",
        dipClusters=gapdir+"/dip_tr_clusters.bed",
        bnBed=expand("calls.bn.{op}.bed",op=shortOps),
        uwbn=expand("bn_supporting_uw.{op}.bed",op=shortOps),
        sup=expand("bn_supporting_uw.{op}.bed",op=shortOps),
        mssmsup=expand("bn_supporting_mssm.{mssm}.{op}.bed",mssm=allMSSMBase.values(),op=shortOps),
        uwBNSupTab="uw_bn_support.tab"        
rule GetBNCalls:
    input:
        bnvcf=config["bnvcf"],        
    output:
        bnBed=expand("calls.bn.{op}.bed",op=shortOps)
    params:
        sd=SD
    shell:"""
zcat {input.bnvcf} | {params.sd}/../sv/utils/variants_vcf_to_bed.py \
  --vcf /dev/stdin --operation deletion --out /dev/stdout --bionano | \
  bioawk -c hdr '{{ if ($svLen >= 50) print;}}' > calls.bn.del.bed

zcat {input.bnvcf} | {params.sd}/../sv/utils/variants_vcf_to_bed.py \
  --vcf /dev/stdin --operation insertion --out /dev/stdout --bionano | \
  bioawk -c hdr '{{ if ($svLen >= 50) print;}}' > calls.bn.ins.bed
"""

rule UWCloseBNCalls:
    input:
        bnBed="calls.bn.{op}.bed",
        uwBed=config["uwsv"]
    output:
        sup="bn_supporting_uw.{op}.bed"
    params:
        longop=lambda wildcards: shortToLong[wildcards.op],
        sd=SNAKEMAKE_DIR,
        maxBNRatio=config["maxRatio"]        
    shell:"""

cat {input.uwBed} | \
  bioawk -c hdr '{{ if ($svType == "{params.longop}") print;}}'  | \
  bedtools intersect -header -loj -a {input.bnBed} -b stdin | \
  bioawk -c hdr '{{ if (NR > 1 && $12 != ".") {{  \
    diff=$svLen-$12; \
     if (diff < 0)  {{ \
       diff=-diff; }} \
     print $1"_"$2"_"$3"\t"$7"_"$8"_"$9"\t"$svLen"\t"diff/$svLen;}} }}' | \
  bioawk -c hdr '{{ if ($4 <= {params.maxBNRatio} ) print }}' > {output.sup}
"""

rule MakeBNSupportTabForUW:
    input:
        uwBed=config["uwsv"],
        bnSup=expand("bn_supporting_uw.{op}.bed",op=shortOps),
    output:
        sup="uw_bn_support.tab"
    params:
        sd=SNAKEMAKE_DIR,
    shell:"""
cat {input.bnSup} | {params.sd}/AddBNSupport.py --svcalls {input.uwBed} --bntable /dev/stdin --out {output.sup}
"""

rule MSSMCloseBNCalls:
    input:
        bnBed="calls.bn.{op}.bed",
        mssmBed="{mssm}.uw.cov",
    output:
        sup="bn_supporting_mssm.{mssm}.{op}.bed"
    params:
        longop=lambda wildcards: shortToLong[wildcards.op],
        sd=SNAKEMAKE_DIR,
        maxBNRatio=config["maxRatio"]        
    shell:"""

cat {input.mssmBed} | \
  bioawk -c hdr '{{ if ($svType == "{params.longop}") print;}}'  | \
  bedtools intersect -header -loj -a {input.bnBed} -b stdin | \
  bioawk -c hdr '{{ if (NR > 1 && $11 != ".") {{  \
    diff=$svLen-$11; \
     if (diff < 0)  {{ \
       diff=-diff; }} \
     print $1"_"$2"_"$3"\t"$7"_"$8"_"$9"\t"$svLen"\t"diff/$svLen;}} }}' | \
  bioawk -c hdr '{{ if ($4 <= {params.maxBNRatio}) print }}' > {output.sup}
"""

rule MSSMCloseBNCallsTable:
    input:
        bnBed=expand("bn_supporting_mssm.{{mssm}}.{op}.bed",op=shortOps),
        mssmBed="{mssm}.uw.cov",        
    output:
        mssmBNTab="{mssm}.sorted.bntab"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
cat {input.bnBed} | {params.sd}/AddBNSupport.py --svcalls {input.mssmBed} --bntable /dev/stdin --out {output.mssmBNTab}
"""
    
       

rule MakeMEIReport:
    input:
        vcf="sv_calls.vcf"
    output:
        mei="mei_report.tab"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
{params.sd}/../sv/PrintMEIReport.py --vcf {input.vcf} --out {output.mei}
"""

################################################################################
##
## BioNano Setup: transfom bionano calls into separate insertion and deletion calls
## with a maximum size.
##
################################################################################

rule MakeBNCalls:
    input:
        bnvcf=config["bnvcf"]
    output:
        bncalls="calls_bn.bed"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
zcat {input.bnvcf} | {params.sd}/../sv/utils/variants_vcf_to_bed.py --vcf /dev/stdin --out {output.bncalls}.tmp --bionano
cat {output.bncalls}.tmp | bioawk -c hdr '{{ if (substr($0,0,1) == "#" || ($svLen >=50 && $svLen < 6000000)) print; }}' > {output.bncalls}
"""
rule MakeSimpleRecirocal:
    input:
        mssm = "mssm.bed.filt",        
        uw   = "uw.bed.filt",
    output:
        simpleReciprocal=expand("reciprocal.{op}.bed",op=shortOps),
    params:
        sge_opts=config["sge_small"],
    shell:"""

bedtools intersect -header -f 0.5 -r \
  -a <( cat {input.mssm} | bioawk -c hdr  '{{ if (NR==1 || $svType == "deletion") print}}') \
  -b <( cat {input.uw} | bioawk -c hdr '{{ if (NR==1 || $svType == "deletion") print}}') -u | wc -l | awk '{{ print $1;}}' > reciprocal.del.bed


bedtools intersect -header -f 0.5 -r \
  -a <( cat {input.mssm} | bioawk -c hdr '{{ if (NR==1 || $svType == "insertion") print}}') \
  -b <( cat {input.uw} | bioawk -c hdr '{{ if (NR==1 || $svType == "insertion") print}}') -u | wc -l | awk '{{ print $1;}}' > reciprocal.ins.bed  
"""

rule MakeWindowRecirocal:
    input:
        mssm = "mssm.bed.filt",        
        uw   = "uw.bed.filt",
    output:
        simpleReciprocal=expand("window.{op}.bed",op=shortOps),
    params:
        sge_opts=config["sge_small"],
    shell:"""

bedtools closest -t first -d \
  -a <( cat {input.mssm} | bioawk -c hdr  '{{ if (NR==1 || $svType == "deletion") print}}') \
  -b <( cat {input.uw} | bioawk -c hdr '{{ if (NR==1 || $svType == "deletion") print}}') | awk '{{ if ($NF < 1000) t+=1 }} END{{ print t"\\t"NR;}}' > window.del.bed


bedtools closest -t first -d \
  -a <( cat {input.mssm} | bioawk -c hdr '{{ if (NR==1 || $svType == "insertion") print}}') \
  -b <( cat {input.uw} | bioawk -c hdr '{{ if (NR==1 || $svType == "insertion") print}}') |  awk '{{ if ($NF < 1000) t+=1 }} END{{ print t"\\t"NR}}' > window.ins.bed
"""


rule SeparateBNCalls:
    input:
        bncalls="calls_bn.bed"
    output:
        bnCallsOp="calls_bn.{op}.bed"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
nf=`head -1 {input.bncalls} | awk '{{ print NF;}}'`
egrep "^#|{wildcards.op}" {input.bncalls} | awk '{{ if (substr($0,0,1) == "#" || $3-$2 < 6000000) print; }}' | bedtools groupby -header -g 1-3 -c 4 -o first -full | cut -f 1-$nf | \
 {params.sd}/../sv/utils/rmdup.py > {output.bnCallsOp}        
"""

shortOpToLongOp={"del": "deletion", "ins" : "insertion"}

################################################################################
##
## MsPAC Setup: transfom bionano calls into separate insertion and deletion calls
## with a maximum size.
##
################################################################################


#
# Copy MSSM calls from source directory
#
rule MakeMSSMLocal:
    input:
        mssm=lambda wildcards: fullMSSM[wildcards.mssm]
    output:
        local="local.{mssm}"
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=01:00:00",
        sn=SNAKEMAKE_DIR,
    shell:"""
bedtools intersect -a {input.mssm} -b {params.sn}/../regions/Regions.Called.bed -u -wa > {output.local}
"""

#
# Transform MSSM calls into UW format so they fit into the QC pipeline
#       
rule MakeMSSMUWpre:
    input:
        mssm="local.{mssm}"
    output:
        mssmUW="{mssm}.uw"
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=01:00:00",
        sd=SNAKEMAKE_DIR
    shell:"""
{params.sd}/../sv/utils/mssm/MSSMToUW.py {input.mssm} --no-reformat --out /dev/stdout | bedtools sort -header | uniq > {output.mssmUW}

"""

#
#  Sort MSSM calls and do some early QC: filter out too small of calls and too large.
#
rule MakeMSSMSort:
    input:
        mssm="{mssm}.uw"
    output:
        mssmSorted="{mssm}.sorted"
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=01:00:00"
    shell:"""
nf=`head -1 {input.mssm} | awk '{{ print NF;}}'`
module unload anaconda; module unload python/3.5.2; module load python/2.7.3; module load bedtools/2.25.0; bedtools sort -header -i {input.mssm} | awk '{{ if (substr($0,0,1) == "#" || ($3-$2 >= 50 && $3-$2 < 500000)) print;}}' | bedtools groupby -g 1-3 -c 4 -o first -full -header | cut -f 1-$nf > {output.mssmSorted}
"""
rule DedupMSSM:
    input:
        sorted="{mssm}.sorted"
    output:
        dedup="{mssm}.sorted.dedup"
    params:
        sd=SNAKEMAKE_DIR,
    shell:"""
cat {input.sorted} | {params.sd}/../sv/utils/rmdup.py > {output.dedup}
"""

rule AnnotateMSSMClusters:
    input:
        sorted="{mssm}.sorted.dedup"
    output:
        mssmClust="{mssm}.sorted.window_clusters"
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=1G -l h_rt=04:00:00 -l disk_free=4G",
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
    shell:"""
cat {input.sorted} | {params.sd}/../sv/utils/ToPoint.sh | \
 bedtools slop -b 250 -g {params.ref}.fai | \
 bedtools sort | \
 bedtools merge -c 2 -o count -i stdin | awk '{{ print $1"\\t"$2"\\t"$3"\\t"$NF;}}' > {output.mssmClust}

"""
        
################################################################################
##
## Add external support for SV calls: BioNano overlap and raw read support
##

#
# First add bionano support for MSSM calls
#

rule AnnotateSourceDataBNOverlap:
    input:
        mssmCalls="{mssm}.sorted",
        bnCalls="calls_bn.bed",
    output:
        mssmBnSupport="bn_support.{mssm}.tab"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
nmssm=`head -1 {input.mssmCalls} | awk '{{ print NF;}}'`
nbm=`head -1 {input.bnCalls} | awk '{{ print NF;}}'`
ratioIndex=$(($nmssm+$nbm+3))

{{ {{ head -1 {input.mssmCalls}; head -1 {input.bnCalls} | {params.sd}/../sv/utils/AddIndex.py _2; }} | paste -s ; bedtools intersect -a {input.mssmCalls} -b {input.bnCalls} -loj; }} | bioawk -c hdr '{{ bnsv=$svLen_2; if (bnsv < 0) {{ bnsv=-1*bnsv; }}; svdiff=bnsv -$svLen; if (svdiff < 0) {{ svdiff =-1*svdiff;}} ; if (substr($0,0,1) == "#") {{ print $0"\\t\\tbnKey\\tbnSize\\tbnOvp"; }} else {{ if ($tStart_2 == "-1" || $svType != $svType_2) {{ print $0"\\tNONE\\t0\\t1";}} else {{ print $0"\\t"$_chrom_2"_"$tStart_2"_"$tEnd_2"\\t"$svLen_2"\\t"svdiff/bnsv; }} }}  }}' | bedtools groupby -header -g 1-3 -c $ratioIndex -o min -full | bioawk -c hdr '{{ print $bnKey"\\t"$bnSize"\\t"$bnOvp;}}' > {output.mssmBnSupport}
"""


##############################################
##
## Add raw read support for MSSM calls
##

##
## Split up file for faster processing.
##
rule SplitMSSM:
    input:
        mssm="{mssm}.sorted.dedup"
    output:
        split=dynamic("{mssm}.sorted.split/gaps.bed.{id}")
    params:
        sge_opts="-pe serial 8 -l mfree=1G -l h_rt=04:00:00",
        sd=SNAKEMAKE_DIR,
        n=config["n"],
        ref=config["ref"]
    shell:"""
module unload anaconda; module unload python/3.5.2; module load python/2.7.3; mkdir -p {wildcards.mssm}.sorted.split; {params.sd}/../sv/utils/SpliceVariantsAndCoverageValidate.py --gaps {input.mssm} --split {params.n} --splitDir {wildcards.mssm}.sorted.split --tmpdir $TMPDIR --ref {params.ref} --blasr {params.sd}/../blasr/alignment/bin/blasr 
"""


##
## Check raw read coverage for MSSM calls
##
rule MakeMSSMSupport:
    input:
        mssm="{mssm}.sorted.split/gaps.bed.{id}",
        reads=config["reads"]
    output:
        cov="{mssm}.sorted.split/gaps_sup.{id}.cov"
    params:
        sge_opts="-pe serial 8 -l mfree=1G -l h_rt=04:00:00",
        ref=config["ref"],
        sd=SNAKEMAKE_DIR
    shell:"""
{params.sd}/../sv/utils/SpliceVariantsAndCoverageValidate.py  --gaps {input.mssm} --reads {input.reads} --out {output.cov} --nproc 1 --ref {params.ref} --tmpdir $TMPDIR --maxSize 50000 --flank 2000 --window 1000 --blasr {params.sd}/../blasr/alignment/bin/blasr
"""

##
## Combine raw read support for MSSM calls.
##
rule CombineMSSMTables:
    input:
        mssmGaps="{mssm}.sorted.split/gaps.bed.{id}",    
        mssmCov="{mssm}.sorted.split/gaps_sup.{id}.cov",
    output:
        mssmComb="{mssm}.sorted.split/gaps_cov.{id}",
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=01:00:00",
        ref=config["ref"]
    shell:
        "paste {input.mssmGaps} {input.mssmCov} > {output.mssmComb}"

rule MergeMSSMSupport:
    input:
        mssmComb=dynamic("{mssm}.sorted.split/gaps_cov.{id}"),    
    output:
        cov="{mssm}.uw.cov"
    params:
        sge_opts="-pe serial 8 -l mfree=1G -l h_rt=01:00:00",
        ref=config["ref"]
    shell:"""
head -1 {wildcards.mssm}.sorted.split/gaps_cov.0 > {output.cov}
grep -hv "^#" {input.mssmComb} | sed '/^\s*$/d' | bedtools sort >> {output.cov}
"""


################################################################################
#
# Now that read coverage and bionano coverage are determined for MSSM calls, filter
# out potentially false positive calls.
#
rule FilterCallsByPBAndBNSupportMSSM:
    input:
        mssmCov="{base}.uw.cov",
        bnSupport="{base}.sorted.bntab"
    output:
        filtMSSM="{base}.uw.bed.filt"
    params:
        sge_opts=config["sge_small"],
        minPbSupport=config["pbSupport"],
        maxBNRatio=config["maxRatio"]
    shell:"""
  paste {input.mssmCov} {input.bnSupport} | \
    bioawk -c hdr \
    '{{ if ( (NR == 1) || ((substr($mssmQC,0,4) == "PASS" || $mssmQC == "hap2_resolved_P" || $mssmQC == "hap1_resolved_P") && \
          (( $tEnd-$tStart < 10000 && $nAlt > {params.minPbSupport}) || $bnstatus == "BN_KEEP" )))  print;}} ' > {output.filtMSSM}
"""

rule MakeMSSMFilteringReport:
    input:
        mssmCov="{base}.uw.cov",
        bnSupport="bn_support.{base}.tab"
    output:
        reportMSSM="{base}.uw.bed.filt.report"
    params:
        sge_opts=config["sge_small"],
        minPbSupport=config["pbSupport"],
        maxBNRatio=config["maxRatio"]
    shell:"""
  paste {input.mssmCov} {input.bnSupport} | \
    bioawk -c hdr \
 'BEGIN {{ qcPass=0; altPass=0;  }} {{ if (NR==1) {{ print "qcPass\\taltPass"; }} \
     else {{ \
       if (((substr($mssmQC,0,4) == "PASS" || $mssmQC == "hap2_resolved_P" || $mssmQC == "hap1_resolved_P" ) )) {{ \
        qcPass = qcPass+1; \
        if (( (( $tEnd-$tStart < 10000 && $nAlt > {params.minPbSupport}) || ($bnKey != "." && $bnOvp < {params.maxBNRatio} )))) {{ \
           altPass = altPass + 1;
        }}\
        }} \
       }} }} \
    END {{ print qcPass"\\t"altPass ; }}' > {output.reportMSSM}
"""

###########################################################
###
### Merging MSPAC haplotypes -- requires merging de novo
### calls and haplotype-specific calls.
### Subtract all haplotype-specific calls from denovos.
### Then add de novos.
###

rule MergeMSSMHapAndDeNovo:
    input:
        mssmh = expand("{mssm}.hap{hap}.bed.uw.bed.filt",mssm=config["sample"],hap=mssmHaps),
        mssmdn = expand("{mssm}.remainingdenovos.bed.uw.bed.filt",mssm=config["sample"])
    output:
        comb = expand("{mssm}.hap{hap}.bed.uw.bed.filt.comb",mssm=config["sample"],hap=mssmHaps),
    params:
        sge_opts=config["sge_small"]
    shell:"""
head -1 {input.mssmh[0]} > {output.comb[0]}
bedtools intersect -v -a {input.mssmdn} -b {input.mssmh[0]} > {input.mssmdn}.no_h1
cat  {input.mssmh[0]} {input.mssmdn}.no_h1  | grep -v "^#" | sort -k1,1 -k2,2n -k3,3n >> {output.comb[0]}
bedtools intersect -v -a {input.mssmdn} -b {input.mssmh[1]} > {input.mssmdn}.no_h2
head -1 {input.mssmh[1]} > {output.comb[1]}
cat   {input.mssmh[1]} {input.mssmdn}.no_h2 | grep -v "^#" | sort -k1,1 -k2,2n -k3,3n >> {output.comb[1]}
"""



################################################################################
##
## Now do tandem repeat cluster filtering of combined calls.  This is quite a few
## steps, and is mirrored from SVQC.
##
## 1. Copy into dir structure that mimics the SVQC dirs.
## 2. Intersect calls iwith tandem_repeats to define clusters.
## 3. Merge tr clusters from haplotypes
## 4. Merge tr clusters with SVQC defined clusters for consistency.
## 5. 
    
rule SetupForFiltering:
    input:
        mssmh = expand("{mssm}.hap{{hap}}.bed.uw.bed.filt.comb",mssm=config["sample"])
    output:
        hapInDir=gapdir+"/hap{hap}/gaps.recalled"
    params:
        sge_opts=config["sge_small"]
    shell:"""
mkdir -p mssm/hap{wildcards.hap}
bedtools sort -header -i {input.mssmh} > {output.hapInDir}
"""
        
rule FindTRClusters:
    input:
        filt=gapdir+"/hap{hap}/gaps.recalled.filt",
    output:
        trClusters=gapdir+"/hap{hap}/tr_clusters.bed",
    params:
        clusterSize=config["tr_cluster_size"],
        sge_opts="-cwd -pe serial 6 -l mfree=1G -l h_rt=12:00:00 -l disk_free=4G",
        contigBed="contigs.h{hap}.fasta.sam.bed",
        sd=SNAKEMAKE_DIR,
        tr=SNAKEMAKE_DIR+"/../regions/tandem_repeats_strs_slop.bed",
    shell:"""
nf=`head -1 {input.filt} | awk '{{ print NF;}}'`
fs=$((nf+1))
fe=$((nf+4))
cat {input.filt} | \
  {params.sd}/../sv/utils/ToPoint.sh | \
  bedtools intersect -header -loj -f 0.9 -a stdin -b {params.tr}  | \
  awk '{{ if (NR==1) {{ print $0"\\ttrChrom\\ttrStart\\ttrEnd\\ttrScore";}} else {{ print; }} }}' | \
  bioawk -c hdr '{{ if (NR > 1 && $trChrom != ".") {{ print $trChrom"\\t"$trStart"\\t"$trEnd;}} }}' |\
  bedtools sort | \
  sort | uniq -c | awk '{{ if ($1 >= {params.clusterSize}) print $2"\\t"$3"\\t"$4"\\t"$1;}}' | \
  bedtools sort > {output.trClusters}
"""

## Step 4, merge from halotypes into diploid
rule MergeTRClusters:
    input:
        trClusters=expand(gapdir+"/hap{hap}/tr_clusters.bed",hap=mssmHaps)
    output:
        dipClusters=gapdir+"/dip_tr_clusters.bed"
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=1G -l h_rt=1:00:00 ",
    shell:"""
cat {input.trClusters} | bedtools sort | bedtools merge > {output.dipClusters}
"""

#
# Remove calls that overlap tandem repeat clusters defined by *either* the MSSM
# tandem clusters, or the UW clusters. This way no calls will overlap the unresolved loci.
#
    
rule SeparateTRClusterSVCalls:
    input:
        filt=gapdir+"/hap{hap}/gaps.recalled.filt",
        trClusters=gapdir+"/dip_tr_clusters.bed",
        uwTrClusters=config["uw_tandem_repeats"]                
    output:
        gapsClust=gapdir+"/hap{hap}/gaps.recalled.clust",
        gapsNoClust=gapdir+"/hap{hap}/gaps.recalled.noclust",
    params:
        sge_opts="-pe serial 1 -l mfree=4G -l h_rt=01:00:00 -l disk_free=4G",
        sd=SNAKEMAKE_DIR,
    shell:"""
cat {input.filt} | {params.sd}/../sv/utils/ToPoint.sh | \
 bedtools intersect -u -header -f 0.9 -a stdin -b {input.trClusters} | \
 {params.sd}/../sv/utils/FromPoint.sh > {output.gapsClust}
            
cat {input.filt} | {params.sd}/../sv/utils/ToPoint.sh | \
 bedtools intersect -v -header -f 0.9 -a stdin -b {input.trClusters} | \
 bedtools intersect -v -header -f 0.9 -a stdin -b {input.uwTrClusters} | \
 {params.sd}/../sv/utils/FromPoint.sh > {output.gapsNoClust}
"""

        

#######################################
# This section is lifted from SVQC.Snakefile, so any changes should be sync'e between here and there.
#
rule AnnotateSVInTR:
    input:
        filt ="mssm/hap{hap}/gaps.recalled.noclust"
    output:
        trNet="mssm/hap{hap}/tr_net.tab"
    params:
        sge_opts="-pe serial 1 -l mfree=4G -l h_rt=01:00:00 -l disk_free=4G",
        tr=SNAKEMAKE_DIR+"/../regions/tandem_repeats_strs_slop.bed",
        sd=SNAKEMAKE_DIR,
    shell:"""
mkdir -p mssm/hap{wildcards.hap}
nf=`head -1 {input.filt} | awk '{{ print NF;}}'`
fs=$(($nf+1))
fe=$(($nf+4))
cat {input.filt} | {params.sd}/../sv/utils/ToPoint.sh | \
  bedtools intersect -header -loj -f 0.9 -a stdin -b {params.tr}  | \
  awk '{{ if (NR==1) {{ print $0"\\ttrChrom\\ttrStart\\ttrEnd\\ttrScore";}} else {{ print; }} }}' | \
  {params.sd}/../sv/utils/ToNet.sh | \
  bioawk -c hdr '{{ if ($trChrom != "." || NR == 1) print;}}' | \
  bedtools groupby  -g $fs-$fe -c 5 -o sum > {output.trNet}
"""

rule FindTRZygosity:
    input:
        trNet=expand("mssm/hap{hap}/tr_net.tab",hap=mssmHaps),
    output:
        trZyg=expand("mssm/hap{hap}/tr_net.tab.zyg",hap=mssmHaps),
        trZygBed=expand("mssm/hap{hap}/tr_net.tab.zyg.bed",hap=mssmHaps),
    params:
        sge_opts="-pe serial 1 -l mfree=4G -l h_rt=01:00:00 -l disk_free=4G",
        tr=SNAKEMAKE_DIR+"/../regions/tandem_repeats_strs_slop.bed",
        sd=SNAKEMAKE_DIR,
    shell:"""
{params.sd}/../sv/utils/AnnotateTRRegionZygosity.py {input.trNet} {output.trZyg}
paste {input.trNet[0]} {output.trZyg[0]} > {output.trZygBed[0]}
paste {input.trNet[1]} {output.trZyg[1]} > {output.trZygBed[1]}
"""

    

rule FilterRecalledGaps:
    input:
        gaps="mssm/hap{hap}/gaps.recalled"
    output:
        filt="mssm/hap{hap}/gaps.recalled.filt"
    params:
        sge_opts="-pe serial 1 -l mfree=4G -l h_rt=01:00:00 -l disk_free=4G",
        inversions=config["inversions"],
    shell:"""
nf=`head -1 {input.gaps} | awk '{{ print NF;}}'`
bedtools intersect -header -f 0.9 -v -a {input.gaps} -b {params.inversions}  | \
bedtools groupby -header -g 1-6 -c 1 -o first -full | cut -f 1-$nf > {output.filt}
"""

    
rule AnnotateTRZygosity:
    input:
        trZygBed="mssm/hap{hap}/tr_net.tab.zyg.bed",
        filt=expand("mssm/hap{{hap}}/gaps.recalled.noclust",mssm=config["sample"]),
    output:
        filtTR="mssm/hap{hap}/gaps.recalled.noclust.tr_bed"
    params:
        sge_opts="-pe serial 1 -l mfree=4G -l h_rt=01:00:00 -l disk_free=4G",
        tr=SNAKEMAKE_DIR+"/../regions/tandem_repeats_strs_slop.bed",
        sd=SNAKEMAKE_DIR,
    shell:"""
nf=`head -1 {input.filt} | awk '{{ print NF;}}'`
fs=$(($nf+1))
fe=$(($nf+4))
cat {input.filt} | {params.sd}/..//sv/utils/ToPoint.sh | \
 bedtools intersect -header -loj -a stdin -b {input.trZygBed} | \
 awk '{{ if (NR==1) {{ print $0"\\ttrChrom\\ttrStart\\ttrEnd\\ttrScore\\ttrExpand\\ttrHap\\ttrHapDiff";}} else {{ print; }}}}' | \
 bedtools groupby -header -g 1-6 -c $fs -o first -full | \
 bioawk -c hdr '{{ print $trHap"\\t"$trChrom"_"$trStart"_"$trEnd"\\t"$trHapDiff;}}' > {output.filtTR}                
"""
    
rule SplitGapsByTR:
    input:
        haps=gapdir+"/hap{hap}/gaps.recalled.noclust",
        filtTR=gapdir+"/hap{hap}/gaps.recalled.noclust.tr_bed",        
    output:
        notr=gapdir+"/hap{hap}/gaps.recalled.noclust.notr",
        homtr=gapdir+"/hap{hap}/gaps.recalled.noclust.homtr",
        hettr=gapdir+"/hap{hap}/gaps.recalled.noclust.hettr",        
    params:
        sge_opts="-cwd -pe serial 8 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        sd=SNAKEMAKE_DIR
    shell:"""
paste {input.haps} {input.filtTR} | bioawk -c hdr '{{ if (NR == 1 || $trHap == ".") print;}}' > {output.notr}
paste {input.haps} {input.filtTR} | bioawk -c hdr '{{ if (NR == 1 || $trHap == "HOM") print;}}' > {output.homtr}
paste {input.haps} {input.filtTR} | bioawk -c hdr '{{ if (NR == 1 || $trHap == "HET") print;}}' > {output.hettr}
"""
    
rule MergeGaps:
    input:
        notr=expand(gapdir+"/hap{hap}/gaps.recalled.noclust.notr",hap=mssmHaps),
        homtr=expand(gapdir+"/hap{hap}/gaps.recalled.noclust.homtr",hap=mssmHaps),
        hettr=expand(gapdir+"/hap{hap}/gaps.recalled.noclust.hettr",hap=mssmHaps),
    output:
        comb="mssm.bed.filt"
    params:
        sge_opts="-cwd -pe serial 8 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        sd=SNAKEMAKE_DIR
    shell:"""
# First combine the calls outside of tandem repeat regions. This is with low threshold for merging.
{params.sd}/../sv/utils/MergeHaplotypesByOperation.sh {input.notr} {output.comb}.not_tr "svType svLen svSeq qName qStart qEnd region nAlt nRef bnKey bnSize bnOvp bnstatus" 0.1

#
# Next combine inside tandem repeat regions that are expected
# to be homozygous. Just one haplotype should be selected here. 
#
cat {input.homtr[0]} | awk '{{ if (NR==1) {{ print $0"\\thap";}} else {{ print $0"\\tHOM";}} }}' | \
 {params.sd}/../sv/utils/Select.py --cols \#chrom tStart tEnd hap svType svLen svSeq qName qStart qEnd region nAlt nRef bnKey bnSize bnOvp bnstatus --out {output.comb}.homtr-0

#
# Finally merge the heterozygous regions, with a moderate threshold on difference.
#
{params.sd}/../sv/utils/MergeHaplotypesByOperation.sh {input.hettr} {output.comb}.hettr "svType svLen svSeq qName qStart qEnd region nAlt nRef bnKey bnSize bnOvp bnstatus" 0.5

# Now combine all calls
head -1 {output.comb}.not_tr > {output.comb}.no_source
cat {output.comb}.not_tr {output.comb}.homtr-0 {output.comb}.hettr | grep -v "^#" | bedtools sort >> {output.comb}.no_source
cat {output.comb}.no_source | awk '{{ if (NR ==1) {{ print $0"\\tsource"; }} else {{ print $0"\\tMSSM";}} }}' > {output.comb}
"""
    
###############################################################
## For both the UW and MSSM calls, 

rule AnnotateBNOverlap:
    input:
        svcalls="{caller}.bed.filt.{op}.bed",
        bncalls="calls_bn.{op}.bed"
    output:
        svBnOvp="bn_overlap.{caller}.bed.filt.{op}.bed",
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
nbn=`head -1 {input.bncalls} | wc -w`
nsv=`head -1 {input.svcalls} | wc -w`
ncol=$(($nbn+$nsv+1))

{{ {{ head -1 {input.bncalls}; head -1 {input.svcalls} | {params.sd}/../sv/utils/AddIndex.py _2; }} | paste -s ; bedtools intersect -a {input.bncalls} -b {input.svcalls} -loj; }} | bioawk -c hdr '{{ bnsv=$svLen; if (bnsv < 0) {{ bnsv=-1*bnsv; }}; svdiff=bnsv -$svLen_2; if (svdiff < 0) {{ svdiff =-1*svdiff;}} ; if (substr($0,0,1) == "#") {{ print $0"\\tbnOvp"; }} else {{ if ($tStart_2 == "-1") {{ print $0"\\t1";}} else {{ print $0"\\t"svdiff/bnsv; }} }}  }}' | bedtools groupby -header -g 1-3 -c $ncol -o min -full | cut -f 1-$ncol > {output.svBnOvp}
"""
    

rule AnnotateBestBNOverlap:
    input:
        uwcalls="bn_overlap.uw.bed.filt.{op}.bed",
        mssmcalls="bn_overlap.mssm.bed.filt.{op}.bed"        
    output:
        optSvCalls="opt_bn_overlap.bed.filt.{op}.bed"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
        maxRatio=config["maxRatio"]
    shell:"""
{params.sd}/CombineBNAnnotatedFiles.py {input.uwcalls} uw {input.mssmcalls} mssm {params.maxRatio} {wildcards.op} > {output.optSvCalls}
"""


rule SelectBestBNOverlap:
    input:
        svCalls="{source}.bed.filt.{op}.bed",
        optCalls="opt_bn_overlap.bed.filt.{op}.bed"
    output:
        svFilt="bn_filt.{source}.{op}.bed"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
{params.sd}/SelectBNSupportedCalls.py --svcalls {input.svCalls} --bntable {input.optCalls} --source {wildcards.source} --out {output.svFilt}
"""

rule SummarizeBestBNOverlap:
    input:
        optCalls="opt_bn_overlap.bed.filt.{op}.bed"
    output:
        optSummary="opt_bn_overlap.bed.filt.{op}.bed.summary"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"]
    shell:"""
mspac=`cat {input.optCalls} | bioawk -c hdr '{{ if ($opSource == "mssm" ) print $bnsvlen; }}' | {params.sd}/stats.py --noround --noheader`
phasedsv=`cat {input.optCalls} | bioawk -c hdr '{{ if ($opSource == "uw" ) print $bnsvlen; }}' | {params.sd}/stats.py --noround --noheader`
echo {params.sample} "mspac" {wildcards.op} $mspac | tr " " "\\t" > {output.optSummary}
echo {params.sample} "phasedsv" {wildcards.op} $phasedsv | tr " " "\\t" >> {output.optSummary}
"""

    
otherSource = { "mssm": "uw", "uw" : "mssm" }



   

rule CombineAnnotatedFixedDeletions:
    input:
        svinsann="fix_del/insertion/insertions.annotated.bed",
        svdelann="fix_del/deletion/deletions.annotated.bed",
    output:
        svfixann="fix_del/sv_calls.ann.bed",
        svfixvcf="fix_del/sv_calls.ann.vcf",
    params:
        sd=SNAKEMAKE_DIR,
        sge_opts=config["sge_small"],
        ref=config["ref"],
    shell:"""
head -1 {input.svdelann} > {output.svfixann}
cat {input.svdelann} {input.svinsann} | bedtools sort >> {output.svfixann}
{params.sd}/../sv/utils/variants_bed_to_vcf.py --bed {output.svfixann} --vcf {output.svfixvcf} --type sv --addci 5 --fields NALT nAlt NREF nRef SVANN svAnn SVREP svRep SVCLASS svClass NTR nTR BN bnKey SOURCE source --source PHASED-SV --info \"##INFO=<ID=NALT,Number=1,Type=Integer,Description=\\\"Number of reads supporting variant\\\">" "##INFO=<ID=NREF,Number=1,Type=Integer,Description=\\\"Number of reads supporting reference\\\">" "##INFO=<ID=SVANN,Number=1,Type=String,Description=\\\"Repeat annotation of variant\\\">" "##INFO=<ID=SVREP,Number=1,Type=Float,Description=\\\"Fraction of SV annotated as mobile element or tandem repeat\\\">"  "##INFO=<ID=SVCLASS,Number=1,Type=String,Description=\\\"General repeat class of variant\\\">"  "##INFO=<ID=NTR,Number=1,Type=Integer,Description=\\\"Number of tandem repeat bases\\\">"  "##INFO=<ID=BN,Number=1,Type=Float,Description=\\\"Overlap with BioNanoGenomics call.\\\">"  "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\\\"Source method of call, with additional information describing haplotype.\\\">" --reference {params.ref}
"""
    
    
rule SummarizeExons:
    input:
        sv="sv_calls.base.bed",
    output:
        exons=config["sample"] + ".exons.{op}.bed",
    params:
        sge_opts=config["sge_small"],
        exonFile=config["exons"]
    shell:"""
geneHeader=`head -1 {params.exonFile}`
egrep "^#|{wildcards.op}" {input.sv} | \
  bioawk -c hdr '{{ if (substr($0,0,1) == "#") \
                   {{ print;}} \
                   else {{ if ($svType == "insertion") {{ $3=$2+1;}} print;}} }}' | \
  tr " " "\\t" | \
  bedtools intersect -header -a stdin -b {params.exonFile} -loj | \
  bedtools groupby -g 1-4 -c 4 -o first -full | \
  awk -v hdr="$header" '{{ if (substr($0,0,1) == "#") print $0"\\t"hdr; }}' > {output.exons}
"""
    

rule MakeMSSMBNAnnotation:
    input:
        mssm="{base}.sorted.dedup",
        bionano=config["bionano"]
    output:
        mssmBN="{base}.uw.bed.bn"
    params:
        sge_opts=config["sge_small"],
        maxRatio=config["maxRatio"]
    shell:
        "~/projects/HGSVG/hgsvg/sv/utils/AnnotateBedWithBioNano.py --bionano {input.bionano} --table {input.mssm} --out {output.mssmBN} --source MSSM --maxRatio {params.maxRatio}"

    
rule MakeUWAnnotation:
    input:
        uw=config["uwsv"],
        bionano=config["bionano"]
    output:
        uwBN="uw.bed.bn"
    params:
        sge_opts=config["sge_small"],
        maxRatio=config["maxRatio"]
    shell:
        "~/projects/HGSVG/hgsvg/sv/utils/AnnotateBedWithBioNano.py --bionano {input.bionano} --table {input.uw} --out {output.uwBN} --source UW --maxRatio {params.maxRatio}"



rule MakeMSSMDistanceReport:
    input:
        mssm="mssm.bed.filt",
        uw="uw.bed.filt"
    output:
        distreport="mssm.bed.filt.distance_report"
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
    shell:"""
nfm=`head -1 {input.mssm} | awk '{{ print NF;}}'`
nfu=`head -1 {input.uw} | awk '{{ print NF;}}'`
nft=$(($nfm+nfu+1))
nelem=`grep -v "^#" {input.mssm} | wc -l`
i=0;
rm -f mssm.shuff;
cut -f 1-3 {input.mssm} > mssm.shuff.bed
nshuff=$(($nfu+4))
while [ $i -lt 10000 ]; do
   bedtools shuffle -i mssm.shuff.bed -g {params.ref}.fai | \
   bedtools sort | \
   bedtools closest -t first -d -a stdin -b {input.uw} | cut -f $nshuff | awk '{{ if ($1 < 1000) print;}}' | wc -l >> mssm.shuff;
   i=$(($i+1))
done
   
nclose=`bedtools closest -t first -d -a {input.mssm} -b {input.uw} | cut -f $nft | awk '{{ if ($1 < 1000) print;}}' | wc -l`
echo -e $nelem"\\t"$nclose > {output.distreport}
"""
    

rule FilterUWCalls:
    input:
        uwsv=config["uwsv"],
        uwbn="uw_bn_support.tab"
    output:
        uwfilt="uw.bed.filt"
    params:
        sge_opts=config["sge_small"],
        pbSupport=config["pbSupport"]
    shell:"""
paste {input.uwsv} {input.uwbn} | bioawk -c hdr '{{ if (NR==1 || $nAlt > {params.pbSupport} || $bnstatus == "BN_KEEP") print;}}' > {output.uwfilt}
"""
            
rule MakeBB:
    input:
        bed="{source}.bed.filt"
    output:
        bb=expand("{{source}}.{op}.bb",op=longOps)
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"]
    shell:"""
module load ucsc
egrep "^#|insertion" {input.bed} | bioawk -c hdr 'BEGIN{{OFS="\\t";}} {{ print $_chrom, $tStart, $tEnd, "{wildcards.source}", 1000, "+", $tStart, $tEnd, "0,0,255";}}' | grep -v "^#" | bedtools sort | {SNAKEMAKE_DIR}/../sv/utils/FixCoordinates.py /dev/stdin /dev/stdout {params.ref}.fai >  {wildcards.source}.insertion.bed
bedToBigBed {wildcards.source}.insertion.bed {params.ref}.fai {wildcards.source}.insertion.bb -type=bed9

egrep "^#|deletion" {input.bed} | bioawk -c hdr 'BEGIN{{colors["HOM"]="255,0,0"; colors["HAP1"]="0,255,0"; colors["HAP0"]="0,0,255"; OFS="\\t";}} {{ print $_chrom, $tStart, $tEnd, "{wildcards.source}", 1000, "+", $tStart, $tEnd, colors[$hap];}}' | grep -v "^#" | bedtools sort | {SNAKEMAKE_DIR}/../sv/utils/FixCoordinates.py /dev/stdin /dev/stdout {params.ref}.fai > {wildcards.source}.deletion.bed
bedToBigBed {wildcards.source}.deletion.bed {params.ref}.fai {wildcards.source}.deletion.bb -type=bed9
"""

rule SeparateByOperation:
    input:
        uw="uw.bed.filt",
        mssm="mssm.bed.filt",
    output:
        svs=expand("{source}.{op}.bed", source=filtered,op=shortOps)
    params:
        sge_opts=config["sge_small"],
    shell:"""
cat {input.uw} | bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $svType == "insertion") print;}}' > {input.uw}.ins.bed
cat {input.uw} | bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $svType == "deletion") print;}}'  > {input.uw}.del.bed
cat {input.mssm} | bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $svType == "insertion") print;}}' > {input.mssm}.ins.bed
cat {input.mssm} | bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $svType == "deletion") print;}}'  > {input.mssm}.del.bed
"""

rule IntersectMSSMAndUW:
    input:
        svs=expand("{source}.{{op}}.bed", source=filtered)
    output:
        isectCounts="{op}.intersect.txt",
    params:
        sge_opts=config["sge_small"],
    shell:"""
rm -f {output.isectCounts}
for frac in 0.0000001 0.1 0.5; do
  nMSSM=`bedtools intersect -a {filtered[0]}.{wildcards.op}.bed -b {filtered[1]}.{wildcards.op}.bed -r -f $frac -wa | wc -l`
  nUW=`bedtools intersect -b {filtered[0]}.{wildcards.op}.bed -a {filtered[1]}.{wildcards.op}.bed -r -f $frac -wa | wc -l`
  echo -e "$frac\\t$nMSSM\\t$nUW" >> {output.isectCounts}
done
"""

rule CountSpliceClusters:
    input:
        sv="{sv}"
    output:
        svClusters="{sv}.cluster-count"
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        bams=config["reads"],
        sd=SNAKEMAKE_DIR
    shell:
        "{params.sd}/../sv/utils/SpliceVariantsAndCoverageValidate.py --gaps {input.sv} --ref {params.ref} --reads {params.bams} --window 250 --flank 1000 --count {output.svClusters} --tmpdir $TMPDIR  --blasr {params.sd}/../blasr/alignment/bin/blasr "


rule SummarizePacBioCovValidation:
    input:
        clusterCount=expand("{uw}.cluster-count", uw=config["uwsv"]) + expand("{mssm}.uw.cluster-count",mssm=allMSSMBase.values()),
        mssmCov=expand("{mssm}.uw.cov",mssm=allMSSMBase.values()),
        uw=config["uwsv"]
    output:
        clusterSummary="cluster-summary.txt",
    params:
        sge_opts=config["sge_small"],
        minPbSupport=config["pbSupport"],
        sample=config["sample"],
        uw=config["uwsv"],
        sd=SNAKEMAKE_DIR        
    shell:"""
nUW=`wc -l {params.uw} | awk '{{ print $1;}}'`
nUWClust=`wc -l {params.uw}.cluster-count | awk '{{ print $1;}}' | awk '{{ print $1;}}'`
nUWPass=`cat {params.uw} | bioawk -c hdr '{{ if ($nAlt > {params.minPbSupport}) print;}}' | wc -l`
echo -e "UW\t$nUW\t$nUWClust\t$nUWPass" > {output.clusterSummary}

nMSSMH1=`wc -l {params.sample}.hap1.bed.uw.bed | awk '{{ print $1;}}'`
nMSSMH1Clust=`wc -l {params.sample}.hap1.bed.uw.cluster-count | awk '{{ print $1;}}'`
nMSSMH1Pass=`cat {params.sample}.hap1.bed.uw.cov | bioawk -c hdr '{{ if ($nAlt > {params.minPbSupport}) print;}}' | wc -l`    
echo -e "MSSMh1\t$nMSSMH1\t$nMSSMH1Clust\t$nMSSMH1Pass" >> {output.clusterSummary}

nMSSMH2=`wc -l {params.sample}.hap2.bed.uw.bed | awk '{{ print $1;}}'`
nMSSMH2Clust=`wc -l {params.sample}.hap2.bed.uw.cluster-count | awk '{{ print $1;}}' `
nMSSMH2Pass=`cat {params.sample}.hap2.bed.uw.cov | bioawk -c hdr '{{ if ($nAlt > {params.minPbSupport}) print;}}' | wc -l`            
echo -e "MSSMh2\t$nMSSMH2\t$nMSSMH2Clust\t$nMSSMH2Pass" >> {output.clusterSummary}

nMSSMDN=`wc -l {params.sample}.remainingdenovos.bed.uw.bed | awk '{{ print $1;}}' `
nMSSMDNClust=`wc -l {params.sample}.remainingdenovos.bed.uw.cluster-count | awk '{{ print $1;}}' `

Rscript {params.sd}/../plotting/plot_inheritance.R --inh {input.merge} --sample {params.sample} --pal {params.pal}
"""

rule SummarizeInheritance:
    input:
        svcalls="sv_calls.annotated.tr_fixed.parents",
    output:
        hethom="het_hom_summary.txt"
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=01:00:00",
        sample=config["sample"],
    shell:"""
nHet=`bioawk -c hdr '{{ if ($hap=="HAP1" || $hap=="HAP0") print; }}' < {input.svcalls} | wc -l`
nHom=`bioawk -c hdr '{{ if ($hap=="HOM") print; }}' < {input.svcalls} | wc -l`
ratio=`echo "scale=3; $nHet/$nHom" | bc`
nHomTP=`cat {input.svcalls} | bioawk -c hdr '{{ if ($hap == "HOM" && $nAlt_fa >= 1 && $nAlt_mo >= 1) print; }}' | wc -l`
nHomShouldBeHet=`cat {input.svcalls}  | bioawk -c hdr '{{ if ($hap == "HOM" && (($nAlt_fa >= 1 && $nAlt_mo < 1) || ($nAlt_fa < 1 && $nAlt_mo >= 1))) print; }}' | wc -l`
nHomFP=`cat {input.svcalls}  | bioawk -c hdr '{{ if ($hap == "HOM" && $nAlt_fa < 1 && $nAlt_mo < 1) print; }}' | wc -l`
nHetTP=`cat {input.svcalls} | bioawk -c hdr '{{ if ( ($hap == "HAP1" || $hap == "HAP0") && ($nAlt_fa >= 1 || $nAlt_mo >=  1)) print; }}' | wc -l`
nHetFP=`cat {input.svcalls} | bioawk -c hdr '{{ if ( ($hap == "HAP1" || $hap == "HAP0") && ($nAlt_fa < 1 && $nAlt_mo < 1))  print; }}' | wc -l`
nNoCall=`bioawk -c hdr '{{ if ($nRef_fa  < 1 && $nAlt_fa < 1 && $nRef_mo < 1 && $nAlt_mo < 1) print; }}' < {input.svcalls} | wc -l`


tnHet=`bioawk -c hdr '{{ if ($is_trf != "TR" && ($hap=="HAP1" || $hap=="HAP0")) print; }}' < {input.svcalls} | wc -l`
tnHom=`bioawk -c hdr '{{ if ($is_trf != "TR" && $hap=="HOM") print; }}' < {input.svcalls} | wc -l`
tratio=`echo "scale=3; $tnHet/$tnHom" | bc`
tnHomTP=`cat {input.svcalls} | bioawk -c hdr '{{ if ($is_trf != "TR" && $hap == "HOM" && $nAlt_fa >= 1 && $nAlt_mo >= 1) print; }}' | wc -l`
tnHomShouldBeHet=`cat {input.svcalls}  | bioawk -c hdr '{{ if ($is_trf ! "TR" && ($hap == "HOM" && (($nAlt_fa >= 1 && $nAlt_mo < 1) || ($nAlt_fa < 1 && $nAlt_mo >= 1)))) print; }}' | wc -l`
tnHomFP=`cat {input.svcalls}  | bioawk -c hdr '{{ if ($is_trf != "TR" && ($hap == "HOM" && $nAlt_fa < 1 && $nAlt_mo < 1)) print; }}' | wc -l`
tnHetTP=`cat {input.svcalls} | bioawk -c hdr '{{ if ( ($is_trf != "TR" && ($hap == "HAP1" || $hap == "HAP0") && ($nAlt_fa >= 1 || $nAlt_mo >=  1))) print; }}' | wc -l`
tnHetFP=`cat {input.svcalls} | bioawk -c hdr '{{ if ($is_trf != "TR" && ( ($hap == "HAP1" || $hap == "HAP0") && ($nAlt_fa < 1 && $nAlt_mo < 1)))  print; }}' | wc -l`
tnNoCall=`bioawk -c hdr '{{ if ($is_trf != "TR" && ($nRef_fa  < 1 && $nAlt_fa < 1 && $nRef_mo < 1 && $nAlt_mo < 1)) print; }}' < {input.svcalls} | wc -l`

echo -e "#condition\\tnHom\\tnHet\\tratio\\thomTP\\thetFPHom\\thomFP\\thetTP\\thetFP\\tNoCall" > {output.hethom}
echo -e "all\\t$nHom\\t$nHet\\t$ratio\\t$nHomTP\\t$nHomShouldBeHet\\t$nHomFP\\t$nHetTP\\t$nHetFP\t$nNoCall" >> {output.hethom}
echo -e "no-trf\\t$tnHom\\t$tnHet\\t$tratio\\t$tnHomTP\\t$tnHomShouldBeHet\\t$tnHomFP\\t$tnHetTP\\t$tnHetFP\t$tnNoCall" >> {output.hethom}        

"""

    
rule AnnotateRecirocalOverlap:
    input:
        uw="uw.bed.filt",
        mssm="mssm.bed.filt",
        svs=expand("{source}.{op}.bed", source=filtered,op=shortOps)        
    output:
        annot=expand("{source}.{op}.overlaps.bed", source=filtered,op=shortOps)
    params:
        sge_opts=config["sge_small"],
    shell:"""

echo "ovp" > {input.uw}.ins.overlaps.bed
{{ {SNAKEMAKE_DIR}/../sv/utils/HeaderMod.py --source {input.uw}.ins.bed {input.mssm}.ins.bed --index; \
bedtools intersect -a {input.uw}.ins.bed -b {input.mssm}.ins.bed -loj -f 0.5 -r  ; }} | \
bedtools groupby -header  -c 4 -o first -full | \
bioawk -c hdr '{{  if ($chrom_2 != ".") {{ print "SHARED";}} else {{ print "UNIQUE";}} }}'| tail -n +2 >> {input.uw}.ins.overlaps.bed

echo "ovp" > {input.uw}.del.overlaps.bed
{{ {SNAKEMAKE_DIR}/../sv/utils/HeaderMod.py --source {input.uw}.del.bed {input.mssm}.del.bed --index; \
bedtools intersect -a {input.uw}.del.bed -b {input.mssm}.del.bed -loj -f 0.5 -r  ; }} | \
bedtools groupby -header  -c 4 -o first -full | \
bioawk -c hdr '{{  if ($chrom_2 != ".") {{ print "SHARED";}} else {{ print "UNIQUE";}} }}'| tail -n +2 >> {input.uw}.del.overlaps.bed


echo "ovp" > {input.mssm}.ins.overlaps.bed
{{ {SNAKEMAKE_DIR}/../sv/utils/HeaderMod.py --source {input.mssm}.ins.bed {input.uw}.ins.bed --index; \
bedtools intersect -a {input.mssm}.ins.bed -b {input.uw}.ins.bed -loj -f 0.5 -r  ; }} | \
bedtools groupby -header  -c 4 -o first -full | \
bioawk -c hdr '{{  if ($chrom_2 != ".") {{ print "SHARED";}} else {{ print "UNIQUE";}} }}'| tail -n +2 >> {input.mssm}.ins.overlaps.bed

echo "ovp" > {input.mssm}.del.overlaps.bed
{{ {SNAKEMAKE_DIR}/../sv/utils/HeaderMod.py --source {input.mssm}.del.bed {input.uw}.del.bed --index; \
bedtools intersect -a {input.mssm}.del.bed -b {input.uw}.del.bed -loj -f 0.5 -r  ; }} | \
bedtools groupby -header  -c 4 -o first -full | \
bioawk -c hdr '{{  if ($chrom_2 != ".") {{ print "SHARED";}} else {{ print "UNIQUE";}} }}'| tail -n +2 >> {input.mssm}.del.overlaps.bed

"""

#
# Rules for merging callsets:
#
#  1. Pick UW or MSSM calls that overlap with BioNano call and less than 20% size difference.
#  2. From remaining calls, if a UW or MSSM call has overlap with bionano, pick that call. Exclude calls from the other dataset that overlap with the bionnano call.
#  3. Pick remaining UW calls, annotated with which ones are 50% reciprocal overlap with MSSM.
#  4. Pick remaining MSSM calls that are under 10kbp, and not closer than 2kbp from a UW call.
#

rule MakeBionano:
    input:
        sv="{source}.bed.filt.{op}.bed"
    output:
        svbn="{source}.bed.filt.{op}.bed.bn"
    params:
        sge_opts=config["sge_small"]
    shell:"""
bioawk -c hdr -t '{{ if ($bnKey != ".") print;}}' < {input.sv} > {output.svbn}
"""

rule AnnotateOtherBioNanoOverlap:
    input:
        sv="{source}.bed.filt.{op}.bed",
        bn="{other}.bed.filt.{op}.bed.bn"
    output:
        svAnnot="{source}.bed.filt.{op}.{other}.bed.filt-bn.bed"
    params:
        sge_opts=config["sge_small"]
    shell:"""
{{ {SNAKEMAKE_DIR}/../sv/utils/HeaderMod.py --source {input.sv} {input.bn} --index ; \
 bedtools intersect -a {input.sv} -b {input.bn} -loj | bedtools groupby -c 4 -o first -full | awk 'BEGIN{{OFS="\\t";}} NF{{NF-=1}};1'; }}  | bioawk -c hdr '{{ if ($chrom_2 != ".") {{ print "{wildcards.other}_BN"; }} else {{ print ".";}} }}' > {output.svAnnot}
"""
rule SplitUWFilt:
    input:
        uw="uw.bed.filt"
    output:
        ops=["uw.bed.filt.ins.filt", "uw.bed.filt.del.filt"]
    params:
        sge_opts=config["sge_small"]
    shell:"""
bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $svType == "insertion") print; }}' < {input.uw} > uw.bed.filt.ins.bed
bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $svType == "deletion") print; }}' < {input.uw} > uw.bed.filt.del.bed
"""

rule AnnotateMSSMDistance:
    input:
        mssm="mssm.bed.filt.{op}.bed",
        uw=expand("uw.bed.filt.{ops}.bed",ops=["del", "ins"]),
    output:
        svAnnot="mssm.bed.filt.{op}.bed.uw-dist",
    params:
        sge_opts=config["sge_small"]
    shell:"""
cat {input.uw[0]} {input.uw[1]} | bedtools sort > {input.uw[0]}.combined
bedtools closest -t first -a {input.mssm} -b {input.uw[0]}.combined -d | \
    awk '{{ if (NR==1) {{ print "distToUW";}} print $NF; }}'> {output.svAnnot}
"""


        
rule AnnotateUWDistance:
    input:
        mssm="mssm.bed.filt.{op}.bed",
        uw="uw.bed.filt.{op}.bed"
    output:
        svAnnot="uw.bed.filt.{op}.bed.mssm-dist",
    params:
        sge_opts=config["sge_small"]
    shell:"""
echo "distToMSSM" > {output.svAnnot}
{{ {SNAKEMAKE_DIR}/../sv/utils/HeaderMod.py --source {input.uw} {input.mssm} --index --append distToMSSM ; \
 bedtools closest -t first -header -a {input.uw} -b {input.mssm} -d | bedtools groupby -c 4 -o first -full; }} | bioawk -c hdr '{{ print $distToMSSM;}}' | tail -n +3 >> {output.svAnnot}
"""

# Merging rules
#
#  1. Pick UW or MSSM calls that overlap with BioNano call and less
#than 20% size difference.
#  2. From remaining calls, if a UW or MSSM call has overlap with
#bionano, pick that call. Exclude calls from the other dataset that
#overlap with the bionnano call.
#  3. Pick remaining UW calls, annotated with which ones are 50%
#reciprocal overlap with MSSM.
#  4. Pick remaining MSSM calls that are under 10kbp, and not closer
#than 2kbp from a UW call.
#  5. Filter out anything  that overlaps an inversion.

    
rule MergeCallsets:
    input:
        uw="uw.bed.filt.{op}.bed",
        mssm="mssm.bed.filt.{op}.bed",
        bnFile=expand("bn_filt.{source}.{{op}}.bed", source=sources),
        uwDist="mssm.bed.filt.{op}.bed.uw-dist",
        mssmDist="uw.bed.filt.{op}.bed.mssm-dist"
    output:
        merged="merged.{op}.bed"
    params:
        sge_opts=config["sge_small"],
        minDist=config["minDist"],
        sd=SD        
    shell:"""

{params.sd}/SelectBestBNKey.py -a {input.uw} -b {input.mssm}
# Remove MSSM calls that have best overlap with bionano from UW. This is the UW callset.
cat {input.mssm}.bnselect | bioawk -c hdr '{{ if ($bnstatus == "BN_KEEP") print $0; }}' > {input.mssm}.keep
    
cat {input.uw}.bnselect | bioawk -c hdr '{{ if ($bnstatus != "BN_OTHER") print;}}' > {input.uw}.no_mssm

nf=`head -1 {input.uw}.no_mssm | awk '{{ print NF;}}'`
bedtools closest -header -d -t first -a {input.uw}.no_mssm -b {input.mssm}.keep | \
  awk '{{ if ($NF > {params.minDist}) print;}}' | \
  cut -f 1-$nf > {input.uw}.far_from_mssm_keep

# Now for the MSSM calls, select either: BN_KEEP, or distToUW > minDist
paste {input.mssm}.bnselect {input.uwDist} | \
  bioawk -c hdr '{{ if ((substr($0,0,1) == "#") || \
    ($distToUW > {params.minDist} || $bnstatus == "BN_KEEP")) print; }}' > {input.mssm}.far_from_uw.bed


{SNAKEMAKE_DIR}/../sv/utils/MergeFiles.py --files  {input.mssm}.far_from_uw.bed {input.uw}.far_from_mssm_keep --out /dev/stdout | bedtools sort -header > {output.merged}
"""


rule CountSourceCalls:
    input:
        uw="uw.bed.filt",
        mssm="mssm.bed.filt",
        bnFile=expand("bn_filt.{source}.{{op}}.bed", source=sources),
        annot=expand("{filt}.{{op}}.overlaps.bed", filt=filtered),
        svbn=expand("{filt}.{{op}}.bed.bn",filt=filtered),
        uwDist="mssm.bed.filt.{op}.bed.uw-dist",
        mssmDist="uw.bed.filt.{op}.bed.mssm-dist"
    output:
        merged="merged.sources.{op}.txt"
    params:
        sge_opts=config["sge_small"],
        minDist=config["minDist"],
        sd=SD
    shell:"""

# Tally calls kept by bionano
mssmBN=`paste mssm.bed.filt.{wildcards.op}.bed bn_filt.mssm.{wildcards.op}.bed | \
   bioawk -c hdr '{{ if ($bnstatus == "BN_KEEP")  print $svLen;}}'  | {params.sd}/../utils/stats.py | tail -n +2 | tr " " "\\t" | cut -f 1-2`
echo -e "MSSM-BN\\t{wildcards.op}\\t$mssmBN"   > {output.merged}


uwBN=`paste uw.bed.filt.{wildcards.op}.bed bn_filt.uw.{wildcards.op}.bed | \
   bioawk -c hdr '{{ if ($bnstatus == "BN_KEEP")  print $svLen;}}'  | {params.sd}/../utils/stats.py | tail -n +2 | tr " " "\\t" | cut -f 1-2`
        
    
echo -e "UW-BN\\t{wildcards.op}\\t$uwBN" >> {output.merged}

echo "mssmmatch" > uw.bed.filt.{wildcards.op}.mssm_match
nf1=`head -1 uw.bed.filt.{wildcards.op}.bed | awk '{{ print NF;}}'`
bedtools intersect  -a  uw.bed.filt.{wildcards.op}.bed -b mssm.bed.filt.{wildcards.op}.bed.keep  -loj | bedtools groupby -g 1-6 -c $nf1 -o first -full | awk '{{ if ($NF == ".") {{ print "NO";}} else {{ print "MATCH";}} }}' >> uw.bed.filt.{wildcards.op}.mssm_match
        
uwNoBn=`cat uw.bed.filt.{wildcards.op}.bed | \
        bioawk -c hdr '{{ if ($bnstatus != "BN_OTHER")  print $svLen;}}'  | {params.sd}/../utils/stats.py | tail -n +2 | tr " " "\\t" | cut -f 1-2`

echo -e "UW-NoBN\\t{wildcards.op}\\t$uwNoBn" >> {output.merged}


mssmRem=`paste mssm.bed.filt.{wildcards.op}.bed mssm.bed.filt.{wildcards.op}.bed.uw-dist  bn_filt.mssm.{wildcards.op}.bed  | bioawk -c hdr '{{ if ($bnstatus == "NA" && $distToUW > {params.minDist}) print $svLen; }}'  | {params.sd}/../utils/stats.py | tail -n +2 | tr " " "\\t" | cut -f 1-2`
echo -e "MSSM-REM\\t{wildcards.op}\\t$mssmRem" >> {output.merged}
"""
    
    
rule CombineMergedToPreInv:
    input:
        merged=expand("merged.{op}.bed", op=shortOps),
    output:
        mergedBed="sv_calls.pre-inv.bed",
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sample=config["sample"],
        sd=SNAKEMAKE_DIR,
    shell:"""
module load numpy/1.11.0    
module load pandas
head -1 merged.ins.bed > {output.mergedBed}.tmp
cat {input.merged} | grep -v "^#" | awk '{{ if ($3-$2 >= 50) print; }}' | bedtools sort >> {output.mergedBed}.tmp
nf=`head -1 {output.mergedBed}.tmp | awk '{{ print NF;}}'`
bedtools groupby -header -g 1-6 -i {output.mergedBed}.tmp -c 1 -o first -full | cut -f 1-$nf > {output.mergedBed}.pre-conflict

{params.sd}/../sv/utils/RemoveDelConflicts.sh  {output.mergedBed}.pre-conflict  {output.mergedBed}

rm -f {output.mergedBed}.tmp
rm -f {output.mergedBed}.pre-conflict

"""

rule ExcludeEventsOverlappingInversions:
    input:
        preInv="sv_calls.pre-inv.bed",
        invHaps="inversions.bed"
    output:
        svCalls="sv_calls.pre-source.bed"
    params:
        sge_opts=config["sge_small"],
    shell:
        "bedtools intersect -header -v -a {input.preInv} -b {input.invHaps} > {output.svCalls}"


rule AddUWSource:
    input:
        svCalls="sv_calls.pre-source.bed",
    output:
        svSource="sv_calls.uw-source.bed"
    params:
        sge_opts=config["sge_small"],
        stitching =[ config["svqc"] + "/hap0/gaps.recalled.sorted",\
                     config["svqc"] + "/hap1/gaps.recalled.sorted"],
        fillin    =[ config["fillin"] + "/hap0/gaps.bed",
                     config["fillin"] + "/hap1/gaps.bed"],
        stitchingFai = [config["stitchingFai"] + "/contigs.h0.fasta.fai", \
                        config["stitchingFai"] + "/contigs.h1.fasta.fai" ],
        fillInFai =[ config["fillInFai"] + "/alignments.h0.bam.fasta.fai", \
                     config["fillInFai"] + "/alignments.h1.bam.fasta.fai" ]
    shell:"""
{SNAKEMAKE_DIR}/../sv/utils/AddOriginalHaplotype.py --sv {input.svCalls} --stitching {params.stitching} --fillin {params.fillin} --out {output.svSource} --stitchingFai {params.stitchingFai} --fillInFai {params.fillInFai} 
"""


rule CountClusters:
    input:
        svUWSource="sv_calls.uw-source.bed",
    output:
        clusteredSVs="cluster-count.tsv",
    params:
        sge_opts=config["sge_small"],
    shell:"""
{SNAKEMAKE_DIR}/../sv/utils/PrintSVClusters.py --gaps {input.svUWSource} --window 500 --out {output.clusteredSVs}
"""

rule AddMSSMSource:
    input:
        svUWSource="sv_calls.uw-source.bed",
        clusteredSVs="cluster-count.tsv",
    output:
        svMSSMSource="sv_calls.base.bed"
    params:
        sge_opts=config["sge_small"],
        mspac=config["mspac"]
    shell:"""
{SNAKEMAKE_DIR}/../sv/utils/AddMsPACSeq.py --sv {input.svUWSource} --mspac {params.mspac} --out {output.svMSSMSource}.tmp
paste {output.svMSSMSource}.tmp {input.clusteredSVs} > {output.svMSSMSource}
rm -f {output.svMSSMSource}.tmp
"""

rule FixDelsInSVCalls:
    input:
        sv="sv_calls.base.bed",
    output:
        svdel="sv_calls.bed.fixed",
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
    shell:"""
{params.sd}/FixDeletionCall.py --svbed {input.sv} --out {output.svdel} --ref {params.ref}
"""

rule AnnotateGapBed:
    input:
        svCalls="sv_calls.bed.fixed"
    output:
        annot="sv_calls.bed.annotated"
    params:
        sge_opts=config["sge_quad"],
        sd=SNAKEMAKE_DIR,
    shell:"""
make -f {params.sd}/../sv/utils/AnnotateGapBed.mak GAPS={input.svCalls}
mv {input.svCalls}.annotated {output.annot}
"""

    
rule AddTRFIntersect:
    input:
        svBase="sv_calls.bed.annotated"
    output:
        trf="in_trf.tab"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
        trf=config["trf_annot"]
    shell:"""
nf1=`head -1 {input.svBase} | awk '{{ print NF;}}'`
nf2=$(($nf1+1))
echo "is_trf" > {output.trf}
cat {input.svBase} | {params.sd}/../sv/utils/ToPoint.sh | bedtools intersect -a stdin -b {params.trf} -f 0.9 -loj | \
bedtools groupby -g 1-6 -c 5 -o first -full | cut -f $nf2 | awk '{{ if ($1 != ".") {{ print "TR"; }} else {{ print ".";}}  }}' >> {output.trf}
"""

rule FixTRAnnotations:
    input:
        annot="sv_calls.bed.annotated",
        trf="in_trf.tab"        
    output:
        fix="sv_calls.bed.annotated.tr_fixed"
    params:
        sge_opts=config["sge_small"],
    shell:"""
paste {input.annot} {input.trf} | bioawk -c hdr '{{ if ($is_trf == "TR" && $svAnn == "NONE") {{ $svAnn = "TandemRepeat";}} print;}}' | tr " " "\t" > {output.fix}        

"""

rule MakeAlt:
    input:
        svCalls="sv_calls.bed.annotated"
    output:
        svCallsAlt="sv_calls.bed.annotated.alt"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
        ref=config["ref"],
        sample=config["sample"]
    shell:"""
{params.sd}/RelabelGenotype.py --bed {input.svCalls} --out {output.svCallsAlt}
"""

rule AddSVQCLoci:
    input:
        svcalls="sv_calls.bed.annotated.tr_fixed",
        loci="../SVQC/tr_clusters.calls.bed"
    output:
        comb="sv_calls.bed.annotated.tr_fixed.loci"
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
    shell:"""
{params.sd}/../sv/utils/MergeFiles.py --files {input.svcalls} {input.loci} --out {output.comb}
cat {output.comb} | bioawk -c hdr '{{ if ($svType == "locus") {{ $svAnn = "locus";}} print;}}' | tr " " "\t" > {output.comb}.relabeled
mv -f {output.comb}.relabeled {output.comb}
"""    
        

rule CombineMergedToVCF:
    input:
        svCalls="sv_calls.bed.annotated.tr_fixed.loci"
    output:
        mergedVCF="sv_calls.vcf",
        mergedVCFgz=config["sample"]+".sv_calls.vcf.gz",
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sample=config["sample"]
    shell:"""

{SNAKEMAKE_DIR}/../sv/utils/variants_bed_to_vcf.py --bed {input.svCalls} --reference {params.ref} --vcf /dev/stdout  --sample {params.sample} --type sv --addci 5 --fields NALT nAlt NREF nRef SVANN svAnn SVREP svRep SVCLASS svClass NTR nTR BN bnKey SOURCE source --source PHASED-SV --info \"##INFO=<ID=NALT,Number=1,Type=Integer,Description=\\\"Number of reads supporting variant\\\">" "##INFO=<ID=NREF,Number=1,Type=Integer,Description=\\\"Number of reads supporting reference\\\">" "##INFO=<ID=SVANN,Number=1,Type=String,Description=\\\"Repeat annotation of variant\\\">" "##INFO=<ID=SVREP,Number=1,Type=Float,Description=\\\"Fraction of SV annotated as mobile element or tandem repeat\\\">"  "##INFO=<ID=SVCLASS,Number=1,Type=String,Description=\\\"General repeat class of variant\\\">"  "##INFO=<ID=NTR,Number=1,Type=Integer,Description=\\\"Number of tandem repeat bases\\\">"  "##INFO=<ID=BN,Number=1,Type=Float,Description=\\\"Overlap with BioNanoGenomics call.\\\">"  "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\\\"Source method of call, with additional information describing haplotype.\\\">" | bedtools sort -header > {output.mergedVCF}

bgzip -c {output.mergedVCF} > {output.mergedVCFgz}
tabix {output.mergedVCFgz}

"""

rule CombineMergedToTracks:
    input:
        svCalls="sv_calls.bed.annotated.tr_fixed.loci"
    output:
        mergedBB=expand("{sample}.{op}.bb",sample=config["sample"],op=shortOps)
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sample=config["sample"]
    shell:"""
cat merged.del.bed | bioawk -c hdr 'BEGIN{{colors["HOM"]="255,0,0"; colors["HAP1"]="0,255,0"; colors["HAP0"]="0,0,255"; OFS="\\t";}} {{ print $_chrom, $tStart, $tEnd, $source, 1000, "+", $tStart, $tEnd, colors[$hap];}}' | grep -v "^#" | bedtools sort | {SNAKEMAKE_DIR}/../sv/utils/FixCoordinates.py /dev/stdin /dev/stdout {params.ref}.fai > {params.sample}.del.bed9
bedToBigBed {params.sample}.del.bed9 {params.ref}.fai {params.sample}.del.bb -type=bed9

cat merged.ins.bed | bioawk -c hdr 'BEGIN{{OFS="\\t";}} {{ print $_chrom, $tStart, $tEnd, $source, 1000, "+", $tStart, $tEnd, "0,0,255";}}' | grep -v "^#" | bedtools sort | {SNAKEMAKE_DIR}/../sv/utils/FixCoordinates.py /dev/stdin /dev/stdout {params.ref}.fai > {params.sample}.ins.bed9

bedToBigBed {params.sample}.ins.bed9 {params.ref}.fai {params.sample}.ins.bb -type=bed9
"""

rule MakeInversions:
    input:
        aln=config["alnDir"] + "/alignments.{hap}.bam"
    output:
        inv=config["alnDir"] + "/inversions/inversions.{hap}.bed"
    params:
        alnDir=config["alnDir"],
        sge_opts="-pe serial 8 -l mfree=3G -l h_rt=24:00:00",
        pbs=config["pbs"],
        ref=config["ref"],
        sample=config["sample"]
    shell:"""

mkdir -p {config.alnDir}
samtools view {input.aln} | {params.pbs}/Inversions/screenInversions /dev/stdin {parmas.ref} {output.inv} -w 10000 -r --noClip -j 12
"""

rule MergeInversions:
    input:
        inv=config["alnDir"] + "/inversions/inversions.{hap}.bed"
    output:
        invMerged="inversions.merged.{hap}.bed"
    params:
        sge_opts=config["sge_small"]
    shell:"""
        echo -e "#chrom\\ttStart\\ttEnd" > {output.invMerged}
        bedtools sort -i {input.inv} | bedtools merge >> {output.invMerged}
"""

rule MergeInversionHaplotypes:
    input:
        invs=expand("inversions.merged.{hap}.bed",hap=haps)
    output:
        invHaps=config["sample"] + ".inversions.UW.bed"
    params:
        sge_opts=config["sge_small"]
    shell:"""
{SNAKEMAKE_DIR}/../sv/utils/MergeHaplotypes.sh {input.invs} {output.invHaps} 
"""



def GetParentBams(w):
    if w == 'fa':
        return config["faBams"]
    else:
        return config["moBams"]
    
rule GenotypeParents:
    input:
        svs="sv_calls.bed.annotated",
        bams=lambda wildcards: GetParentBams(wildcards.pa)
    output:
        parent="parent.cov.{pa}.bed"
    params:
        sge_opts="-pe serial 8 -l mfree=1G -l h_rt=08:00:00",
        ref=config["ref"],
        faBam=config["faBams"],
        moBam=config["moBams"],
        sd=SNAKEMAKE_DIR
    shell:"""
for i in $(seq 0 7); do
  {SNAKEMAKE_DIR}/../sv/utils/SpliceVariantsAndCoverageValidate.py --gaps {input.svs} --reads {input.bams} --out {output.parent}.$i --nproc 1 --ref {params.ref} --separate-haplotypes  --blasr {params.sd}/../blasr/alignment/bin/blasr  --paIndex $i --paNumber 8 &
done
wait
head -1 {output.parent}.0 > {output.parent}
for i in $(seq 0 7); do
tail -n +2 {output.parent}.$i >> {output.parent}
rm -f {output.parent}.i
done
                
perl -pi -e "s/#region/region_{wildcards.pa}/g" {output.parent}
perl -pi -e "s/nAlt/nAlt_{wildcards.pa}/g" {output.parent}
perl -pi -e "s/nRef/nRef_{wildcards.pa}/g" {output.parent}
    
"""

rule MergeInheritance:
    input:
        parents=expand("parent.cov.{pa}.bed",pa=parents),        
        sv="sv_calls.bed.annotated.tr_fixed",
    output:
        merge="sv_calls.annotated.tr_fixed.parents",
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=01:00:00",
    shell:"paste {input.sv} {input.parents}  > {output.merge}"

rule RenderGenotypes:
    input:
        merge="sv_calls.annotated.tr_fixed.parents",        
    output:
        plot=config["sample"]+".parental_coverage.pdf",
        plotnotrf=config["sample"]+".parental_coverage.no_trf.pdf"        
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=01:00:00",
        sample=config["sample"],
        pal=config["pal"],
        sd=SNAKEMAKE_DIR
    shell:"""

Rscript {params.sd}/../plotting/plot_inheritance.R --inh {input.merge} --sample {params.sample} --pal {params.pal}
"""


rule MakeBNCompare:
    input:
        bnvcf=config["bnvcf"],
        svcalls="sv_calls.annotated.tr_fixed.parents",
    output:
        bnCompare=["sv_calls.bed.ins.bn.bed", "sv_calls.bed.del.bn.bed", "sv_calls.bed.by-pb.bn.bed"],
    params:
        sge_opts="-pe serial 8 -l mfree=1G -l h_rt=08:00:00",
    shell:"""
{SNAKEMAKE_DIR}/../sv/utils/IntersectBioNanoWithUW.sh {input.bnvcf} {input.svcalls}

mv sv_calls.annotated.tr_fixed.parents.ins.bn.bed sv_calls.bed.ins.bn.bed
mv sv_calls.annotated.tr_fixed.parents.del.bn.bed sv_calls.bed.del.bn.bed
mv sv_calls.annotated.tr_fixed.parents.by-pb.bn.bed sv_calls.bed.by-pb.bn.bed
"""

rule MakeBNPlots:
    input:
        bnCompare=["sv_calls.bed.ins.bn.bed", "sv_calls.bed.ins.bn.bed", "sv_calls.bed.by-pb.bn.bed"]
    output:
        bnPlot=expand("BioNano.PacBio.Agreement.{sample}.pdf",sample=config["sample"])
    params:
        sge_opts="-pe serial 8 -l mfree=1G -l h_rt=08:00:00",
        sample=config["sample"]
    shell:"""
module load R/latest; Rscript {SNAKEMAKE_DIR}/../plotting/plot_bionano.R --pbsv sv_calls.bed.by-pb.bn.bed --bndel sv_calls.bed.del.bn.bed --bnins sv_calls.bed.ins.bn.bed --sample {params.sample}
"""
