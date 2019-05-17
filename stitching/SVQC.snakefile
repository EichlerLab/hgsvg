config = {}
config["haps"] = ["0", "1"]
config["recall_bin"]= 100

configfile: "phasedsv.json"


haps=config["haps"]
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
import json
grid_file=open("grid.json")
gridOpts= json.load(grid_file)
for key in gridOpts:
    config[key] = gridOpts[key]
    

ops=["insertion", "deletion"]
dirs=["SVQC", "fill-in"]
shell.prefix("source ./config.sh; ")
gapdir="SVQC"
config["ngmlr_cutoff"] = 1000

if "contigs" in config:
    contigs=config["contigs"]
else:
    contigs = { "0" : "contigs.h0.fasta", "1" : "contigs.h1.fasta" }

if "do_recall" not in config:
    config["do_recall"] = "yes"



    
rule all:
    input:
        stitching=expand("stitching_hap_gaps/hap{hap}/gaps.bed", hap=haps),
        sortedBed=expand("SVQC/hap{hap}/gaps.sorted",hap=haps),
        orig_clusters=expand("SVQC/hap{hap}/gaps.sorted.window_clusters",hap=haps),
        bb=expand("SVQC/hap{hap}/{op}s.bb",hap=haps,op=["insertion","deletion"]),
        mergedRetained=expand("SVQC/hap{hap}/indels.svqc.bed",hap=haps),
        filtMerged="merged/filt_sv_calls.bed",
        filtVCF="merged/filt_sv_calls.vcf",
        mergedBed="merged/sv_calls.bed",
        mergedVCF="merged/sv_calls.vcf",
	sampleVCF="merged/" + config["sample"] + ".sv_calls.vcf.gz",
        mergedRetainedVCF=expand("SVQC/hap{hap}/indels.svqc.raw.vcf",hap=haps),
        normVCF=expand("SVQC/hap{hap}/indels.svqc.norm.vcf",hap=haps),
        normBED=expand("SVQC/hap{hap}/indels.svqc.norm.bed",hap=haps),
        normOpBed=expand("SVQC/hap{hap}/indels.svqc.norm.{op}.bed",hap=haps,op=ops),
        normOpBedSup=expand("SVQC/hap{hap}/indels.svqc.norm.{op}.bed.support",hap=haps,op=ops),
        dipIndels="SVQC/diploid/indels.bed",
        dipOpIndels=expand("SVQC/diploid/indels.{op}.bed",op=ops),
        dipVCF=expand("SVQC/diploid/{sample}.indels.vcf",sample=config["sample"]),
        localIndel=expand("fill-in/hap{hap}/indels.{op}.bed",hap=haps,op=ops),
        fillInAll=expand("fill-in/hap{hap}/gaps.all_regions.bed",hap=haps),
        fillInHap=expand("fill-in/hap{hap}/gaps.bed",hap=haps),
        fillInCov=expand("fill-in/hap{hap}/gaps.bed.cov",hap=haps),
        fillInFilt=expand("fill-in/hap{hap}/gaps.bed.support",hap=haps),
        fillInRegion=expand("fill-in/hap{hap}/region.bed",hap=haps),
        fillInFiltNoClust=expand("fill-in/hap{hap}/gaps.bed.support.noclust",hap=haps),        
        fillInFiltClusters=expand("fill-in/hap{hap}/tr_clusters.bed",hap=haps),        
        fillInDiploid="fill-in/diploid/sv_calls.bed.support", 
        fillInComb=expand("fill-in/hap{hap}/comb.bed",hap=haps),
        fillInFiltVCCF=expand("fill-in/hap{hap}/"+config["sample"]+".hap{hap}.vcf",hap=haps),
        fillInFiltVCCFGZ=expand("fill-in/hap{hap}/"+config["sample"]+".hap{hap}.vcf.gz",hap=haps),
        dipTrClusters=gapdir+"/dip_tr_clusters.bed",
        dipAllTrClusters=gapdir+"/dip_tr_clusters.bed.all",        
        reCalled=expand("{gapdir}/hap{hap}/gaps.recalled",gapdir=gapdir,hap=haps),
        reCalledSorted=expand("{gapdir}/hap{hap}/gaps.recalled.sorted",gapdir=gapdir,hap=haps),
        reCalledFiltered=expand("{gapdir}/hap{hap}/gaps.recalled.filt",gapdir=gapdir,hap=haps),
        reCalledFilteredClusters=expand("{gapdir}/hap{hap}/gaps.recalled.filt.clusters",gapdir=gapdir,hap=haps),
        indelBed=expand("{gapdir}/hap{hap}/indels.recalled.bed",gapdir=gapdir,hap=haps),
        filtSVCalls=expand("{gapdir}/diploid/sv_calls.bed.filt",gapdir=gapdir),
        svCalls=expand("{gapdir}/diploid/sv_calls.bed",gapdir=gapdir),
        retainedIndels=expand("{gapdir}/hap{hap}/indels.retained.bed",gapdir=gapdir, hap=haps),
        recalledRgions=expand("{gapdir}/hap{hap}/regions.recalled.ref",gapdir=gapdir, hap=haps),
        recalledSortedRegions=expand("{gapdir}/hap{hap}/regions.recalled.ref.sorted",gapdir=gapdir, hap=haps),
        trNet=expand("{gapdir}/hap{hap}/tr_net.tab",gapdir=gapdir,hap=haps),
        trZyg=expand("{gapdir}/hap{hap}/tr_net.tab.zyg",gapdir=gapdir,hap=haps),
        trZygBed=expand("{gapdir}/hap{hap}/tr_net.tab.zyg.bed",gapdir=gapdir,hap=haps),   
        filtTR=expand("{gapdir}/hap{hap}/gaps.recalled.noclust.tr_bed",gapdir=gapdir,hap=haps),
        notr=expand("{gapdir}/hap{hap}/gaps.recalled.noclust.notr",gapdir=gapdir,hap=haps),
        hettr=expand("{gapdir}/hap{hap}/gaps.recalled.noclust.hettr",gapdir=gapdir,hap=haps),
        trClusters=expand("{gapdir}/hap{hap}/tr_clusters.bed",gapdir=gapdir,hap=haps),
        allTrClusters=expand("{gapdir}/hap{hap}/tr_clusters.bed.all",gapdir=gapdir,hap=haps),        
        trClusterLifted=expand("{gapdir}/hap{hap}/tr_clusters.bed.to_asm",gapdir=gapdir,hap=haps),
        trClusterFasta=expand("{gapdir}/hap{hap}/tr_clusters.bed.fasta",gapdir=gapdir,hap=haps),
        gapsClust=expand("{gapdir}/hap{hap}/gaps.recalled.clust",gapdir=gapdir,hap=haps),
        gapsNoClust=expand("{gapdir}/hap{hap}/gaps.recalled.noclust",gapdir=gapdir,hap=haps),
        trClusterCalls=gapdir+"/tr_clusters.calls.bed",
        indels=expand("SVQC/hap{hap}/indels_local.bed",hap=haps),
        indels_op=expand("SVQC/hap{hap}/indels_local.{op}.bed",hap=haps, op=ops),
        fillInSVCoverage=expand("fill-in/hap{hap}/gaps.bed.readcov",hap=haps),
        hapSVCoverage=expand("SVQC/hap{hap}/gaps.recalled.sorted.readcov",hap=haps)
rule HapSVCov:
    input:
        calls="SVQC/hap{hap}/gaps.recalled.sorted"
    output:
        svCoverage="SVQC/hap{hap}/gaps.recalled.sorted.readcov"
    params:
        grid_opts=config["grid_large"],
        sd=SNAKEMAKE_DIR,
        fofn=config["bams"],
        depth=config["depth"]
    shell:"""
{params.sd}/../sv/utils/SVCoverage.py --fofn {params.fofn} --calls {input.calls} --out {output.svCoverage} --nproc 1  --op DEL  --header --window 1000 --default {params.depth}
"""

rule FillInHapSVCov:
    input:
        calls="fill-in/hap{hap}/gaps.bed"
    output:
        svCoverage="fill-in/hap{hap}/gaps.bed.readcov"
    params:
        grid_opts=config["grid_large"],
        sd=SNAKEMAKE_DIR,
        fofn=config["bams"],
        depth=config["depth"]
    shell:"""
{params.sd}/../sv/utils/SVCoverage.py --fofn {params.fofn} --calls {input.calls} --out {output.svCoverage} --nproc 1  --op DEL --header --window 1000 --default {params.depth}
"""


rule CombinedHaps:
    input:
        fillInGaps="fill-in/hap{hap}/gaps.bed.support",
        asmGaps=gapdir+"/hap{hap}/gaps.recalled.filt",
    output:
        comb="fill-in/hap{hap}/comb.bed",
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR,
        ref=config["ref"],
    shell:"""
bedtools intersect -header -v -a {input.fillInGaps} -b {input.asmGaps} > {input.fillInGaps}.no_asm
{params.sd}/../sv/utils/MergeFiles.py --files {input.asmGaps} {input.fillInGaps}.no_asm --out /dev/stdout | \
  bedtools sort -header > {output.comb}
"""

rule FilterCombinedHaps:
    input:
        comb="fill-in/hap{hap}/comb.bed",
    output:
        filt="fill-in/hap{hap}/comb.bed.filt",
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR,
        ref=config["ref"],
        depth=int(float(config["depth"])*0.60)
    shell:"""
cat {input.comb} | bioawk -c hdr '{{ if (NR==1 || ($nAlt > 3 || ($svType == "deletion" && $coverage < {params.depth}))) {{ print;}} }}' > {output.filt}
"""

rule CombinedHapsToVCF:
    input:
        filt="fill-in/hap{hap}/comb.bed.filt",
    output:
        vcf="fill-in/hap{hap}/"+config["sample"]+".hap{hap}.vcf"
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR,
        ref=config["ref"],
        depth=int(float(config["depth"])*0.60),
        sample=config["sample"]
    shell:"""
{SNAKEMAKE_DIR}/../sv/utils/variants_bed_to_vcf.py --bed {input.filt} --ref {params.ref} --sample {params.sample} --type sv --seq --vcf /dev/stdout --fields NALT nAlt NREF nRef | bedtools sort -header > {output.vcf}

"""

rule IndexHapVCF:
    input:
        vcf="fill-in/hap{hap}/"+config["sample"]+".hap{hap}.vcf"
    output:
        gz="fill-in/hap{hap}/"+config["sample"]+".hap{hap}.vcf.gz",
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
bgzip -c {input.vcf} > {output.gz}
tabix {output.gz}
"""

rule CreateClusterCalls:
    input:
        reCalledFiltered=gapdir+"/hap{hap}/gaps.recalled.filt",
        contigBed=lambda wildcards: contigs[wildcards.hap]+".sam.bed"

    output:
        reCalledFilteredClusters=gapdir+"/hap{hap}/gaps.recalled.filt.clusters"
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR,
        ref=config["ref"],
        clusterSize=config["tr_cluster_size"]
    shell:"""
cat {input.reCalledFiltered} | {params.sd}/../sv/utils/ToPoint.sh | \
   bedtools slop -g {params.ref}.fai -b 1000 -i stdin |
   bedtools sort | \
   bedtools merge > {output.reCalledFilteredClusters}.regions

bedtools intersect -a {output.reCalledFilteredClusters}.regions -b {input.reCalledFiltered} -wa | \
  bedtools groupby -g 1-3 -c 2 -o count | \
  awk '{{ if ($NF >= {params.clusterSize}) print;}}' > {output.reCalledFilteredClusters}

rm -f {output.reCalledFilteredClusters}.regions  
"""

rule CreateTRClusterCalls:
    input:
        regions=gapdir+"/dip_tr_clusters.bed",
        lifted=expand(gapdir+"/hap{hap}/tr_clusters.bed.to_asm",hap=haps),
        fasta=expand(gapdir+"/hap{hap}/tr_clusters.bed.fasta",hap=haps),
    output:
        calls=gapdir+"/tr_clusters.calls.bed"
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR
    shell:"""
{params.sd}/../sv/utils/CombineVariableLoci.py  --regions {input.regions} --lifted {input.lifted}  --sequences {input.fasta} --outFile {output.calls}
"""
    
rule MergeAllTRClusters:
    input:
        clusters=expand(gapdir+"/hap{hap}/tr_clusters.bed.all", hap=haps),
    output:
        dip=gapdir+"/dip_tr_clusters.bed.all",
    params:
        grid_opts=config["grid_small"],
    shell:"""
cat {input.clusters} | bedtools sort | bedtools merge > {output.dip}
"""
    
rule MergeTRClusters:
    input:
        clusters=expand(gapdir+"/hap{hap}/tr_clusters.bed", hap=haps),
        allClusters=expand(gapdir+"/hap{hap}/gaps.recalled.filt.clusters",hap=haps)        
    output:
        dip=gapdir+"/dip_tr_clusters.bed",
    params:
        grid_opts=config["grid_small"],
    shell:"""
cat {input.clusters} {input.allClusters} | bedtools sort | bedtools merge > {output.dip}
"""
    
rule SplitGaps:
    input:
        gaps=gapdir+"/hap{hap}/gaps.sorted",
        asm=lambda wildcards: contigs[wildcards.hap]
    output:
        splitGaps=dynamic(gapdir+"/hap{hap}/split/gaps.bed.{id}")
    params:
        n=config["recall_bin"],
        grid_opts=config["grid_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
        gd=gapdir
    shell:"""
mkdir -p {params.gd}/hap{wildcards.hap}/split;
{params.sd}/RecallRegionsInGapBed.py --asm {input.asm} \
                                     --ref {params.ref} \
                                     --gaps {input.gaps} \
                                     --split {params.n} \
                                     --splitDir {params.gd}/hap{wildcards.hap}/split

"""


    
rule RecallGaps:
    input:
        gaps=gapdir+"/hap{hap}/split/gaps.bed.{id}",
        asm=lambda wildcards: contigs[wildcards.hap]
    output:
        recalled=gapdir+"/hap{hap}/split/gaps.recalled.{id}",
        refRegions=gapdir+"/hap{hap}/split/regions.recalled.ref.{id}",
    params:
        grid_opts=config["grid_manycore"],
        ref=config["ref"],
        ngmlr_cutoff=config["ngmlr_cutoff"],
        do_recall=config["do_recall"],
        sd=SNAKEMAKE_DIR,
        gapdir=gapdir
    shell:"""


mkdir -p {params.gapdir}/hap{wildcards.hap};
mkdir -p {params.gapdir}/hap{wildcards.hap}/indels;
if [ {params.do_recall} == "yes" ]; then
  {params.sd}/RecallRegionsInGapBed.py --asm {input.asm} --ref {params.ref} --gaps {input.gaps} --out {output.recalled} --nproc 2 --refRegions {params.gapdir}/hap{wildcards.hap}/split/regions.recalled.ref.{wildcards.id} --ngmlr {params.ngmlr_cutoff} --indels {output.recalled}.indel.bed --indelDir {params.gapdir}/hap{wildcards.hap}/indels --blasr {params.sd}/../blasr/alignment/bin/blasr
else
  cp {input.gaps} {output.recalled}
  touch {output.refRegions}
fi

"""

rule MergeRecalledIndels:
    input:
        gapBed=gapdir+"/hap{hap}/gaps.recalled"
    output:
        indelBed=gapdir+"/hap{hap}/indels.recalled.bed"
    params:
        grid_opts=config["grid_small"],
        gapdir=gapdir
    shell:"""
cat {params.gapdir}/hap{wildcards.hap}/indels/* | head -1 > {output.indelBed}
nf=`head {output.indelBed} | awk '{{ print NF;}}'`
cat {params.gapdir}/hap{wildcards.hap}/indels/* | grep -v "^#" | \
cut -f 1-$nf | \
awk -v fields=$nf '{{if (NF==fields) print;}}' | bedtools sort >> {output.indelBed}
"""

rule MergeRetainedAndRecalledIndels:
    input:
        indelBed=gapdir+"/hap{hap}/indels.recalled.bed",
        retainedIndels=gapdir+"/hap{hap}/indels.retained.bed",
    output:
        mergedRetained=gapdir+"/hap{hap}/indels.svqc.bed"
    params:
        grid_opts=config["grid_small"],
    shell:"""
head -1 {input.indelBed} | cut -f 1-6 > {output.mergedRetained}
cat {input.indelBed} {input.retainedIndels} | cut -f 1-6 | grep -v "^#" | bedtools sort >> {output.mergedRetained}
"""
    
rule SplicedPBSupport:
    input:
        gaps=gapdir+"/hap{hap}/split/gaps.recalled.{id}"
    params:
        grid_opts=config["grid_manycore"],
        ref=config["ref"],
        bams=config["bams"],
        sd=SNAKEMAKE_DIR,
        gapdir=gapdir
    output:
        pbSupport=gapdir+"/hap{hap}/split/gaps.recalled_support.{id}"
    shell:"""

bedtools sort -header -i {input.gaps} > {input.gaps}.tmp
mv -f {input.gaps}.tmp {input.gaps}
mkdir -p {params.gapdir}/hap{wildcards.hap};
{params.sd}/../sv/utils/SpliceVariantsAndCoverageValidate.py --gaps {input.gaps} --ref {params.ref} --reads {params.bams} --window 1000 --flank 2000 --out {output.pbSupport} --nproc 1 --blasr {params.sd}/../blasr/alignment/bin/blasr 
"""        

rule AddSplicedPBSupport:
    input:
        sup=gapdir+"/hap{hap}/split/gaps.recalled_support.{id}",
        gaps=gapdir+"/hap{hap}/split/gaps.recalled.{id}"        
    params:
        grid_opts=config["grid_manycore"],
        ref=config["ref"],
        bams=config["bams"]
    output:
        pbSupport=gapdir+"/hap{hap}/split/gaps.gaps_recalled_support.{id}"
    shell:
        "paste  {input.gaps} {input.sup} > {output.pbSupport}"

rule MergeRecallGaps:
    input:
        gaps=dynamic(gapdir+"/hap{hap}/split/gaps.gaps_recalled_support.{id}"),
        refRegions=dynamic(gapdir+"/hap{hap}/split/regions.recalled.ref.{id}")
    output:
        recalled=gapdir+"/hap{hap}/gaps.recalled",
        recalledRegions=gapdir+"/hap{hap}/regions.recalled.ref"        
    params:
        grid_opts=config["grid_manycore"],
        gapdir=gapdir
    shell:"""
# Create the header, but some of the pasted headers have extra #'s in them, so
# remove all but the first
head -1 {params.gapdir}/hap{wildcards.hap}/split/gaps.gaps_recalled_support.0 | \
    tr -d "#" | awk '{{ print "#"$0; }}'> {output.recalled}
cat {input.gaps} | grep -v "^#" | bedtools sort >> {output.recalled}
cat {input.refRegions} > {output.recalledRegions}
"""



rule SortRecallGaps:
    input:
        recalled="{gapdir}/hap{hap}/gaps.recalled"
    output:
        recalledSorted="{gapdir}/hap{hap}/gaps.recalled.sorted"
    params:
        grid_opts=config["grid_small"],
    shell:"""
bedtools sort -header -i {input.recalled} > {output.recalledSorted}
"""


rule FilterRecalledGaps:
    input:
        gaps="{gapdir}/hap{hap}/gaps.recalled.sorted",
        readcov="{gapdir}/hap{hap}/gaps.recalled.sorted.readcov",
    output:
        filt="{gapdir}/hap{hap}/gaps.recalled.filt"
    params:
        grid_opts=config["grid_small"],
        inversions=config["inversions"],
    shell:"""
nf=`paste {input.gaps} {input.readcov} | head -1 | awk '{{ print NF;}}'`
paste {input.gaps} {input.readcov} | bedtools intersect -header -f 0.9 -v -a stdin -b {params.inversions}  | \
bedtools groupby -header -g 1-5 -c 1 -o first -full | cut -f 1-$nf > {output.filt}
"""

rule AnnotateAllTRClusters:
    input:
        filt="{gapdir}/hap{hap}/gaps.recalled.filt",
    output:
        trClusters="{gapdir}/hap{hap}/tr_clusters.bed.all",
    params:
        clusterSize=config["tr_cluster_size"],
        grid_opts=config["grid_manycore"],
        contigBed=lambda wildcards: contigs[wildcards.hap] + ".sam.bed",
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
  sort | uniq -c | awk '{{ if ($1 >= 2 ) print $2"\\t"$3"\\t"$4"\\t"$1;}}' | \
  bedtools sort > {output.trClusters}
"""

    
rule FindTRClusters:
    input:
        filt="{gapdir}/hap{hap}/gaps.recalled.filt",
    output:
        trClusters="{gapdir}/hap{hap}/tr_clusters.bed",
    params:
        clusterSize=config["tr_cluster_size"],
        grid_opts=config["grid_manycore"],
        contigBed=lambda wildcards: contigs[wildcards.hap] + ".sam.bed",
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


        
rule SeparateTRClusterSVCalls:
    input:
        filt=gapdir+"/hap{hap}/gaps.recalled.filt",
        trClusters=gapdir+"/dip_tr_clusters.bed",
    output:
        gapsClust=gapdir+"/hap{hap}/gaps.recalled.clust",
        gapsNoClust=gapdir+"/hap{hap}/gaps.recalled.noclust",
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
cat {input.filt} | {params.sd}/../sv/utils/ToPoint.sh | \
 bedtools intersect -u -header -f 0.9 -a stdin -b {input.trClusters} | \
 {params.sd}/../sv/utils/FromPoint.sh | tr " " "\t" > {output.gapsClust}
            
cat {input.filt} | {params.sd}/../sv/utils/ToPoint.sh | \
 bedtools intersect -v -header -f 0.9 -a stdin -b {input.trClusters} | \
 {params.sd}/../sv/utils/FromPoint.sh | tr " " "\t" > {output.gapsNoClust}
"""        
    
rule AnnotateSVInTR:
    input:
        filt="{gapdir}/hap{hap}/gaps.recalled.noclust"
    output:
        trNet="{gapdir}/hap{hap}/tr_net.tab"
    params:
        grid_opts=config["grid_small"],
        tr=SNAKEMAKE_DIR+"/../regions/tandem_repeats_strs_slop.bed",
        sd=SNAKEMAKE_DIR,
    shell:"""

nf=`head -1 {input.filt} | awk '{{ print NF;}}'`
fs=$((nf+1))
fe=$((nf+4))
cat {input.filt} | {params.sd}/../sv/utils/ToPoint.sh | \
  bedtools intersect -header -loj -f 0.9 -a stdin -b {params.tr}  | \
  awk '{{ if (NR==1) {{ print $0"\\ttrChrom\\ttrStart\\ttrEnd\\ttrScore";}} else {{ print; }} }}' | \
  {params.sd}/../sv/utils/ToNet.sh | \
  bioawk -c hdr '{{ if ($trChrom != "." || NR == 1) print;}}' | \
  bedtools groupby  -g $fs-$fe -c 5 -o sum > {output.trNet}
"""

    
rule FindTRZygosity:
    input:
        trNet=expand("{gapdir}/hap{hap}/tr_net.tab",gapdir=gapdir,hap=haps),
    output:
        trZyg=expand("{gapdir}/hap{hap}/tr_net.tab.zyg",gapdir=gapdir,hap=haps),
        trZygBed=expand("{gapdir}/hap{hap}/tr_net.tab.zyg.bed",gapdir=gapdir,hap=haps),        
    params:
        grid_opts=config["grid_small"],
        tr=SNAKEMAKE_DIR+"/../regions/tandem_repeats_strs_slop",
        sd=SNAKEMAKE_DIR,
    shell:"""
{params.sd}/../sv/utils/AnnotateTRRegionZygosity.py {input.trNet} {output.trZyg}
paste {input.trNet[0]} {output.trZyg[0]} > {output.trZygBed[0]}
paste {input.trNet[1]} {output.trZyg[1]} > {output.trZygBed[1]}
"""

rule AnnotateTRZygosity:
    input:
        trZygBed=gapdir+"/hap{hap}/tr_net.tab.zyg.bed",
        filt=gapdir+"/hap{hap}/gaps.recalled.noclust"
    output:
        filtTR=gapdir+"/hap{hap}/gaps.recalled.noclust.tr_bed"
    params:
        grid_opts=config["grid_small"],
        tr=SNAKEMAKE_DIR+"/../regions/tandem_repeats_strs_slop.bed",
        sd=SNAKEMAKE_DIR,
    shell:"""
nf=`head -1 {input.filt} | awk '{{ print NF;}}'`
fc=$(($nf+1))
cat {input.filt} | {params.sd}/../sv/utils/ToPoint.sh | \
 bedtools intersect -header -loj -a stdin -b {input.trZygBed} | \
 awk '{{ if (NR==1) {{ print $0"\\ttrChrom\\ttrStart\\ttrEnd\\ttrScore\\ttrExpand\\ttrHap\\ttrHapDiff";}} else {{ print; }}}}' | \
 bedtools groupby -header -g 1-5 -c $fc -o first -full | \
 bioawk -c hdr '{{ print $trHap"\\t"$trHapDiff;}}' > {output.filtTR}
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
        grid_opts=config["grid_manycore"],
        sd=SNAKEMAKE_DIR
    shell:"""
paste {input.haps} {input.filtTR} | bioawk -c hdr '{{ if (NR == 1 || $trHap == ".") print;}}' > {output.notr}
paste {input.haps} {input.filtTR} | bioawk -c hdr '{{ if (NR == 1 || $trHap == "HOM") print;}}' > {output.homtr}
paste {input.haps} {input.filtTR} | bioawk -c hdr '{{ if (NR == 1 || $trHap == "HET") print;}}' > {output.hettr}
"""

rule FilterGaps:
    input:
        comb=gapdir+"/diploid/sv_calls.bed",
    output:
        filt=gapdir+"/diploid/sv_calls.bed.filt",        
    params:
        grid_opts=config["grid_small"],
        cov=int(config["depth"])*0.6
    shell:"""
cat {input.comb} | bioawk -c hdr '{{ if (NR==1 || $nAlt > 3 || ($svType == "deletion" && $coverage < {params.cov})) {{ print;}} }}' > {output.filt}
"""
 

rule MergeGaps:
    input:
        notr=expand(gapdir+"/hap{hap}/gaps.recalled.noclust.notr", hap=haps),
        homtr=expand(gapdir+"/hap{hap}/gaps.recalled.noclust.homtr", hap=haps),
        hettr=expand(gapdir+"/hap{hap}/gaps.recalled.noclust.hettr", hap=haps),
        trClusters=gapdir+"/dip_tr_clusters.bed",
    output:
        comb=gapdir+"/diploid/sv_calls.bed"
    params:
        grid_opts=config["grid_manycore"],
        sd=SNAKEMAKE_DIR
    shell:"""
# First combine the calls outside of tandem repeat regions. This is with low threshold for merging.
{params.sd}/../sv/utils/MergeHaplotypesByOperation.sh {input.notr} {output.comb}.not_tr "svType svLen svSeq qName qStart qEnd region nAlt nRef coverage" 0.1

#
# Next combine inside tandem repeat regions that are expected
# to be homozygous. Just one haplotype should be selected here. 
#
cat {input.homtr[0]} | awk '{{ if (NR==1) {{ print $0"\\thap";}} else {{ print $0"\\tHOM";}} }}' | tr " " "\t" | \
 {params.sd}/../sv/utils/Select.py --cols \#chrom tStart tEnd hap svType svLen svSeq qName qStart qEnd region nAlt nRef coverage --out {output.comb}.homtr-0

#
# Finally merge the heterozygous regions, with a moderate threshold on difference.
#
{params.sd}/../sv/utils/MergeHaplotypesByOperation.sh {input.hettr} {output.comb}.hettr "svType svLen svSeq qName qStart qEnd region nAlt nRef coverage" 0.5

# Now combine all calls
head -1 {output.comb}.not_tr > {output.comb}.pre-filter
cat {output.comb}.not_tr {output.comb}.homtr-0 {output.comb}.hettr | grep -v "^#" | bedtools sort >> {output.comb}.pre-filter

{params.sd}/../sv/utils/RemoveDelConflicts.sh {output.comb}.pre-filter {output.comb}

#rm -f {output.comb}.pre-filter
        
"""


rule MakeBB:
    input:
        gaps="SVQC/hap{hap}/gaps.recalled.noclust"
    output:
        bbs=expand("SVQC/hap{{hap}}/{op}s.bb",op=["deletion","insertion"]),
    params:
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
        grid_opts=config["grid_small"]
    shell:"""
cat {input.gaps} | bioawk -c hdr '{{ if (NR==1 || ($svType == "deletion" && $nAlt > 3) ) print;}}' | {params.sd}/../sv/utils/FixCoordinates.py /dev/stdin {input}.del {params.ref}.fai
cat {input.gaps} | bioawk -c hdr '{{ if (NR==1 || ($svType == "insertion" && $nAlt > 3) ) print;}}' | {params.sd}/../sv/utils/FixCoordinates.py /dev/stdin {input}.ins {params.ref}.fai
{params.sd}/../sv/utils/GapBedToBed6.py {input.gaps}.del {input.gaps}.del6
{params.sd}/../sv/utils/GapBedToBed6.py {input.gaps}.ins {input.gaps}.ins6
bedToBigBed {input.gaps}.del6 {params.ref}.fai {output.bbs[0]} -type=bed6
bedToBigBed {input.gaps}.ins6 {params.ref}.fai {output.bbs[1]} -type=bed6
"""

rule SortRecalledRegions:
    input:
        recalled="SVQC/hap{hap}/regions.recalled.ref",
    output:
        rsorted="SVQC/hap{hap}/regions.recalled.ref.sorted",
    params:
        grid_opts=config["grid_small"]
    shell:
        "bedtools sort -i {input.recalled} > {output.rsorted}"

        
    
rule CollectRetainedIndels:
    input:
        stitchIndels="stitching_hap_gaps/hap{hap}/indels.norm.bed",
        recalledRegions="SVQC/hap{hap}/regions.recalled.ref.sorted"
    output:
        retainedIndels="SVQC/hap{hap}/indels.retained.bed",
    params:
        grid_opts=config["grid_small"],
    shell:
        "bedtools intersect -v -a {input.stitchIndels} -b {input.recalledRegions} > {output.retainedIndels} "


    
rule SortGaps:
    input:
        stitching="stitching_hap_gaps/hap{hap}/gaps.bed"
    output:
        svqc="SVQC/hap{hap}/gaps.sorted"
    params:
        grid_opts=config["grid_small"],
    shell:"""
mkdir -p SVQC/hap{wildcards.hap}
nf=`head -1 {input.stitching} | awk '{{ print NF;}}'`
bedtools sort -header -i {input.stitching} | \
 bedtools groupby -header -g 1-5 -c 1 -o first -full | cut -f 1-$nf > {output.svqc}
"""

rule CountGapClusters:
    input:
        svqc="SVQC/hap{hap}/gaps.sorted"
    output:
        orig_clusters="SVQC/hap{hap}/gaps.sorted.window_clusters"
    params:
        grid_opts=config["grid_small"],
        ref=config["ref"],        
        sd=SNAKEMAKE_DIR,        
    shell:"""
cat {input.svqc} | {params.sd}/../sv/utils/ToPoint.sh | \
 bedtools slop -b 250 -g {params.ref}.fai | \
 bedtools sort | \
 bedtools merge -c 2 -o count -i stdin | awk '{{ print $1"\\t"$2"\\t"$3"\\t"$NF;}}' > {output.orig_clusters}
"""
rule ConvertIndelBedToVCF:
    input:
        mergedRetained="SVQC/hap{hap}/indels.svqc.bed"
    output:
        mergedRetainedVCF="SVQC/hap{hap}/indels.svqc.raw.vcf"
    params:
        grid_opts=config["grid_small"],
        ref=config["ref"],
        sample=config["sample"]
    shell:"""
{SNAKEMAKE_DIR}/../sv/utils/variants_bed_to_vcf.py --bed {input.mergedRetained} --ref {params.ref} --sample {params.sample} --type indel --vcf /dev/stdout | bedtools sort -header > {output.mergedRetainedVCF}
"""

rule NormIndelVCF:
    input:
        rawVCF="SVQC/hap{hap}/indels.svqc.raw.vcf"
    output:
        normVCF="SVQC/hap{hap}/indels.svqc.norm.vcf"
    params:
        grid_opts=config["grid_small"],
        ref=config["ref"],
    shell:"""
vt normalize -r {params.ref} -o {output.normVCF} {input.rawVCF}
"""

rule NormIndelVCFToBed:
    input:
        vcf="SVQC/hap{hap}/indels.svqc.norm.vcf"
    output:
        bed="SVQC/hap{hap}/indels.svqc.norm.bed"
    params:
        grid_opts=config["grid_small"],
        ref=config["ref"],
    shell:"""
{SNAKEMAKE_DIR}/../sv/utils/variants_vcf_to_bed.py --vcf {input.vcf} --out {output.bed}
"""
    
rule SplitNormBed:
    input:
        allbed="SVQC/hap{hap}/indels.svqc.norm.bed"
    output:
        opbed="SVQC/hap{hap}/indels.svqc.norm.{op}.bed"
    params:
        grid_opts=config["grid_small"],
        ref=config["ref"],
    shell:"""
egrep "^#|{wildcards.op}" {input.allbed} > {output.opbed}
"""

rule LocalAssemblyBasedIndels:
    input:
        aln="alignments.h{hap}.bam"
    output:
        indels="SVQC/hap{hap}/indels_local.bed"
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR,
        ref=config["ref"],
    shell:"""
samtools view {input.aln}| {params.sd}/../sv/utils/PrintGaps.py {params.ref} /dev/stdin --maxLength 50 --minLength 2 --outFile /dev/stdout  |\
  bedtools sort -header > {output.indels}
"""

rule SplitSupport:
    input:
        indels="SVQC/hap{hap}/indels_local.bed"
    output:
        indels_op=expand("SVQC/hap{{hap}}/indels_local.{op}.bed",op=ops)
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR,
        ref=config["ref"],
    shell:"""
for op in insertion deletion; do 
  cat {input.indels} | bioawk -c hdr -v sv=$op '{{ if (NR ==1 || $svType == sv) print ;}}' > SVQC/hap{wildcards.hap}/indels_local.$op.bed
done
"""
    
rule AddSupport:
    input:
        opbed="SVQC/hap{hap}/indels.svqc.norm.{op}.bed",
        locbed="SVQC/hap{hap}/indels_local.{op}.bed"
    output:
        opsupport="SVQC/hap{hap}/indels.svqc.norm.{op}.bed.support",
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR,
        ref=config["ref"],
    shell:"""
cat {input.opbed} | \
 bioawk -c hdr '{{ if (NR==1) {{ print "#oChrom\\toStart\\toEnd\\t"$0;}} else {{ print $1"\\t"$2"\\t"$3"\\t"$0}}}}' |  \
 bedtools slop -header -i stdin -g {params.ref}.fai  -b 200 | \
 bedtools intersect -header -a stdin -b {input.locbed} -loj  |\
 awk '{{ if ($NF != ".") print;}}' |\
 awk '{{ if (NR == 1) {{ print $0"\\tqChrom\\tqStart\\tqEnd\\tqop\\tqsvlen\\tqsvseq\\tqsvtsd\\tqsvasm"; }} else {{ print $0;}}}}' | \
  {params.sd}/../indels/CountLocalIndelSupport.py | bioawk -c hdr '{{ if (NR == 1 || $locsup > 1) print;}}' > {output.opsupport}
"""

rule MergeSupport:
    input:
       opsupport=expand("SVQC/hap{{hap}}/indels.svqc.norm.{op}.bed.support",op=ops)
    output:
       mergedopsupport="SVQC/hap{hap}/indels.svqc.norm.bed.support",
    params:
        grid_opts=config["grid_small"],
    shell:"""
head -1 {input.opsupport[0]} > {output.mergedopsupport}
cat {input.opsupport} | grep -v "^#" | bedtools sort >> {output.mergedopsupport}
"""

rule HapIndelBedToDiploidBedOp:
    input:
        indels=expand("SVQC/hap{hap}/indels.svqc.norm.{{op}}.bed.support",hap=haps)
    output:
        dipOpIndels="SVQC/diploid/indels.{op}.bed"
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR
    shell:"""
{params.sd}/../indels/MergeHaplotypes.sh {input.indels} {output.dipOpIndels}
"""


rule HapIndelBedToDiploidBed:
    input:
        dipopindels=expand("SVQC/diploid/indels.{op}.bed",op=ops)
    output:
        dipIndels="SVQC/diploid/indels.bed"
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR
    shell:"""
head -1 {input.dipopindels[0]} > {output.dipIndels}
cat {input.dipopindels} | grep -v "^#" | bedtools sort >> {output.dipIndels}
"""

rule IndelBedToVCF:
    input:
        bed="SVQC/diploid/indels.bed"
    output:
        vcf=expand("SVQC/diploid/{sample}.indels.vcf",sample=config["sample"])
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR,
        ref=config["ref"],
        sample=config["sample"]
    shell:"""

{SNAKEMAKE_DIR}/../sv/utils/variants_bed_to_vcf.py --bed {input.bed} --ref {params.ref} --sample {params.sample} --type indel --vcf /dev/stdout | bedtools sort -header > {output.vcf}
        bgzip -c {output.vcf} > {output.vcf}.gz
        tabix {output.vcf}.gz
"""
    
    
rule LiftTRClusters:
    input:
        trClusters=gapdir+"/dip_tr_clusters.bed",
    output:
        trClusterLifted="SVQC/hap{hap}/tr_clusters.bed.to_asm",
    params:
        clusterSize=config["tr_cluster_size"],
        grid_opts=config["grid_manycore"],
        contigSam=lambda wildcards: contigs[wildcards.hap] + ".sam",
        sd=SNAKEMAKE_DIR
    shell:"""
{params.sd}/../mcutils/src/samLiftover {params.contigSam} {input.trClusters} {output.trClusterLifted} --dir 1 --printNA
"""
    
rule MakeTRClusterFasta:
    input:
        trClusterLifted="SVQC/hap{hap}/tr_clusters.bed.to_asm",
    output:
        trClusterFasta="SVQC/hap{hap}/tr_clusters.bed.fasta",
    params:
        clusterSize=config["tr_cluster_size"],
        grid_opts=config["grid_manycore"],
        contigs=lambda wildcards: contigs[wildcards.hap]
    shell:"""
ctg=`head -1 {params.contigs} | cut -f 1 | tr -d ">"`
samtools faidx {params.contigs} `cat {input.trClusterLifted} | awk -v ctg=$ctg '{{ if ($1 != "NA") {{ print $1":"$2"-"$3;}} else {{ print ctg":1-1";}} }}' | tr "\\n" " "` > {output.trClusterFasta}  
"""
    
    
rule MakeSVOpBins:
    input:
        indelbed="SVQC/bins.bed"
    output:
        separateBins=expand("SVQC/bins.{op}.bed", op=ops)
    params:
        grid_opts=config["grid_manycore"],
        ref=config["ref"],
        prindel="/net/eichler/vol5/home/mchaisso/projects/mcst/prindel",
        covWalks=SNAKEMAKE_DIR+"/../sv/utils/CovBinsWalks.py"
    shell:"""
{params.covWalks} SVQC/bins.bed --op ins --out SVQC/bins.insertion.bed;
{params.covWalks} SVQC/bins.bed --op del --out SVQC/bins.deletion.bed;
"""

    
#rule AnnotateSVGaps:
#    input:
#         combined=expand("{{dir}}/hap{hap}/gaps.recalled.sorted",hap=haps)
#    output:
#         calls="{dir}/diploid/sv_calls.bed"
#    params:
#        grid_opts="-cwd -pe serial 8 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
#        sd=SNAKEMAKE_DIR
#    shell:"""
#
#make -f {params.sd}/SVQC.DiploidAnnotation.mak GAPS=gaps.bed.support DIR={wildcards.dir}
#"""
#

rule FillInRegions:
    input:
        contigBed=lambda wildcards: contigs[wildcards.hap] + ".sam.bed",        
    output:
        fillInRegion="fill-in/hap{hap}/region.bed"
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR
    shell:"""
bedtools subtract -a {params.sd}/../regions/Regions.Called.bed -b {input.contigBed} > {output.fillInRegion}
"""
    

rule MakeFillIn:
    input:
        contigBam="alignments.h{hap}.bam",
        fillInRegion="fill-in/hap{hap}/region.bed"
    output:
        fillInHap="fill-in/hap{hap}/gaps.all.bed"
    params:
        grid_opts=config["grid_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR
#
# Grab SVs from fill-in. Since it's possible the alignments came from different
# regions from what are allowed, filter out from the called regions.
#
    shell:"""
mkdir -p fill-in/hap{wildcards.hap}
samtools view {input.contigBam} `cat {input.fillInRegion} | awk '{{ print $1":"$2"-"$3;}}' | tr "\\n" " "` | \
sort -k1,1 -k2,2n | \
uniq | \
{params.sd}/../sv/utils/PrintGaps.py {params.ref} /dev/stdin --maxMasked 10 --minAlignmentLength 30000 --minContigLength 30000 --condense 20  | \
  bedtools sort -header | \
  {params.sd}/../sv/utils/rmdup.py > {output.fillInHap}
"""


rule KeepNotCoveredFillIn:
    input:
        fillInHap="fill-in/hap{hap}/gaps.all.bed",
        contigBed=lambda wildcards: contigs[wildcards.hap] + ".sam.bed"        
    output:
        fillIn="fill-in/hap{hap}/gaps.all_regions.bed"
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR
    shell:"""
bedtools intersect -header -v -a {input.fillInHap} -b {input.contigBed} > {output.fillIn}
"""

rule RemoveHeterochromatic:
    input:
        fillInAll="fill-in/hap{hap}/gaps.all_regions.bed"
    output:
        fillIn="fill-in/hap{hap}/gaps.bed"
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR
    shell:"""
bedtools intersect -header -a {input.fillInAll} -b {params.sd}/../regions/Regions.Called.bed -wa -u > {output.fillIn}
"""

rule FillInIndel:
    input:
        contigBam="alignments.h{hap}.bam",
        fillInRegion="fill-in/hap{hap}/region.bed"        
    output:
        indel="fill-in/hap{hap}/indels.{op}.bed",
    params:
        grid_opts=config["grid_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR        
    shell:"""
  mkdir -p fill-in/hap{wildcards.hap}
  samtools view {input.contigBam} `cat {input.fillInRegion} | awk '{{ print $1":"$2"-"$3;}}' | tr "\n" " "` | \
  sort -k1,1 -k2,2n | \
  uniq | \
    {params.sd}/../sv/utils/PrintGaps.py {params.ref} /dev/stdin --minLength 2 --maxLength 50 | \
  grep {wildcards.op} | bedtools sort -header | \
  bedtools intersect -header -a stdin -b {params.sd}/../regions/Regions.Called.bed -wa -u > {output.indel}
"""
        
        
rule FillInSupport:
    input:
        fillIn="fill-in/hap{hap}/gaps.bed"
    output:
        fillInCov="fill-in/hap{hap}/gaps.bed.cov"
    params:
        grid_opts=config["grid_small"],
        ref=config["ref"],
        bams=config["bams"],
        sd=SNAKEMAKE_DIR
    shell:"""
{params.sd}/../sv/utils/SpliceVariantsAndCoverageValidate.py --blasr {params.sd}/../blasr/alignment/bin/blasr --gaps {input.fillIn} --ref {params.ref} --reads {params.bams} --window 250 --flank 1000 --out {output.fillInCov} --nproc 1 --blasr {params.sd}/../blasr/alignment/bin/blasr 
"""



rule FilteredFillIn:
    input:
        fillInCov="fill-in/hap{hap}/gaps.bed.cov",
        fillInReadCov="fill-in/hap{hap}/gaps.bed.readcov",
        fillIn="fill-in/hap{hap}/gaps.bed"
    output:
        fillInFilt="fill-in/hap{hap}/gaps.bed.support"
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR,
        cov=int(config["depth"])*0.6
    shell:"""
paste {input.fillIn} {input.fillInCov} {input.fillInReadCov} | bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $nAlt > 3 || ($svType == "deletion" && $coverage < {params.cov}) ) print;}}' > {output.fillInFilt}
"""

rule FillInTRClusters:
    input:
        fillInFilt="fill-in/hap{hap}/gaps.bed.support"
    output:
        fillInFiltClusters="fill-in/hap{hap}/tr_clusters.bed"
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR,
        cluster_count=config["tr_cluster_size"]
    shell:"""
cat {input.fillInFilt} | \
{params.sd}/../sv/utils/ToPoint.sh | \
bedtools intersect -a stdin -b {params.sd}/../regions/tandem_repeats_strs_slop.bed -wb | \
awk '{{ a=NF; print $(a-3)"\\t"$(a-2)"\\t"$(a-1)"\\t"$a;}}' | \
bedtools sort | \
bedtools groupby -g 1-3 -c 2 -o count -full | \
awk -vcc={params.cluster_count} '{{ if ($NF >= cc) print;}}' > {output.fillInFiltClusters}
"""

rule RemoveClusteredSVsFromFilt:
    input:
        fillInFilt="fill-in/hap{hap}/gaps.bed.support",
        fillInFiltClusters="fill-in/hap{hap}/tr_clusters.bed"
    output:
        fillInFilt="fill-in/hap{hap}/gaps.bed.support.noclust",
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
cat {input.fillInFilt} | \
{params.sd}/../sv/utils/ToPoint.sh | \
bedtools intersect -v -header -a stdin -b {input.fillInFiltClusters} | \
{params.sd}/../sv/utils/FromPoint.sh > {output.fillInFilt}
"""

rule FilteredDiploid:
    input:
        fillInFilt=expand("fill-in/hap{hap}/gaps.bed.support.noclust",hap=haps)
    output:
        fillInDiploid="fill-in/diploid/sv_calls.bed.support"
    params:
        grid_opts=config["grid_small"],
        sd=SNAKEMAKE_DIR
    shell:"""
{params.sd}/../sv/utils/MergeHaplotypesByOperation.sh {input.fillInFilt} {output.fillInDiploid}.pre "svType svLen svSeq qName qStart qEnd region nAlt nRef coverage"

{params.sd}/../sv/utils/RemoveDelConflicts.sh {output.fillInDiploid}.pre {output.fillInDiploid}
rm -f {output.fillInDiploid}.pre

"""
    
rule MakeFiltMergedBed:
    input:
        fillInFilt="fill-in/diploid/sv_calls.bed.support",
        filtSVCalls="SVQC/diploid/sv_calls.bed.filt",
    output:
        filtMerged="merged/filt_sv_calls.bed",
    params:
        grid_opts=config["grid_small"],
    shell:"""
mkdir -p merged
# Make sure fill in does not overlap with stitching.
bedtools intersect -header -v -a {input.fillInFilt} -b {input.filtSVCalls} > merged/fill-in.bed.filt

# Prepare name
head -1 merged/fill-in.bed.filt | awk '{{ print $0"\\tsource";}}' > merged/fill-in.bed.filt.src
tail -n +2 merged/fill-in.bed.filt | awk '{{ print $0"\\tlocal";}}' >>  merged/fill-in.bed.filt.src

head -1 {input.filtSVCalls} | awk '{{ print $0"\\tsource";}}' > merged/svqc.bed.src.filt
tail -n +2 {input.filtSVCalls} |  awk '{{ print $0"\\tstitching";}}' >> merged/svqc.bed.src.filt

 {SNAKEMAKE_DIR}/../sv/utils/Select.py --table merged/fill-in.bed.filt.src  --out merged/fill-in.bed.filt.subset --cols `head -1 {input.filtSVCalls}` source
 {SNAKEMAKE_DIR}/../sv/utils/MergeFiles.py --files merged/svqc.bed.src.filt merged/fill-in.bed.filt.src | bedtools sort -header > {output.filtMerged}
    
"""
 



rule MakeMergedBed:
    input:
        svqcBed="SVQC/diploid/sv_calls.bed",
        fillinBed="fill-in/diploid/sv_calls.bed.support"
    output:
        mergedBed="merged/sv_calls.bed"
    params:
        grid_opts=config["grid_small"],
    shell:"""
mkdir -p merged
bedtools intersect -header -v -a {input.fillinBed} -b {input.svqcBed} > merged/fill-in.bed
head -1 merged/fill-in.bed | awk '{{ print $0"\\tsource";}}' > merged/fill-in.bed.src
tail -n +2 merged/fill-in.bed | awk '{{ print $0"\\tlocal";}}' >>  merged/fill-in.bed.src

head -1 {input.svqcBed} | awk '{{ print $0"\\tsource";}}' > merged/svqc.bed.src 
tail -n +2 {input.svqcBed} |  awk '{{ print $0"\\tstitching";}}' >> merged/svqc.bed.src

 {SNAKEMAKE_DIR}/../sv/utils/Select.py --table merged/fill-in.bed.src  --out merged/fill-in.bed.subset --cols `head -1 {input.svqcBed}` source
 {SNAKEMAKE_DIR}/../sv/utils/MergeFiles.py --files merged/svqc.bed.src merged/fill-in.bed.src | bedtools sort -header > {output.mergedBed}
    
"""
#    .format(SNAKEMAKE_DIR, SNAKEMAKE_DIR)

    
rule MakeMergedFiltVCF:
    input:
        mergedBed="merged/filt_sv_calls.bed",
    output:
        mergedVCF="merged/filt_sv_calls.vcf",
    params:
        grid_opts=config["grid_small"],
        ref=config["ref"], 
        sample=config["sample"]
    shell:"""
{SNAKEMAKE_DIR}/../sv/utils/variants_bed_to_vcf.py --bed {input.mergedBed} --ref {params.ref} --sample {params.sample} --type sv --seq --vcf /dev/stdout --fields NALT nAlt NREF nRef SRC source | bedtools sort -header > {output.mergedVCF}"""
    
rule MakeMergedVCF:
    input:
        mergedBed="merged/sv_calls.bed"
    output:
        mergedVCF="merged/sv_calls.vcf"
    params:
        grid_opts=config["grid_small"],
        ref=config["ref"], 
        sample=config["sample"]
    shell:"""
{SNAKEMAKE_DIR}/../sv/utils/variants_bed_to_vcf.py --bed {input.mergedBed} --ref {params.ref} --sample {params.sample} --type sv --vcf /dev/stdout --fields NALT nAlt NREF nRef SRC source | bedtools sort -header > {output.mergedVCF}"""


rule MakeSampleVCF:
    input:
        mergedVCF="merged/sv_calls.vcf"
    output:
        sampleVCF="merged/"+config["sample"]+".sv_calls.vcf.gz"
    params:
        grid_opts=config["grid_small"],
    shell:
        "bgzip -c {input.mergedVCF} > {output.sampleVCF}; tabix {output.sampleVCF}"
        
    
