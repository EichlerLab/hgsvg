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

configfile: "sv_support.json"


SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
SD=SNAKEMAKE_DIR
shell.prefix(". {SNAKEMAKE_DIR}/config.sh; ")

cwd=os.getcwd()

svTypes=["DEL", "INS"]
opToOperation={"del": "deletion", "ins": "insertion"}
pbSVTypeMap={"DEL": "deletion", "INS": "insertion"}
bnSVTypeMap={"DEL": "del", "INS": "ins"}

bioNanoOps = { "DEL": "bncalls.DEL.bed", "INS": "bncalls.INS.bed" }

haps=["0", "1"]

annotationTables=  [SD+"/SegDups/hg38.conservedElts.v2.bed",
                    SD+"/SegDups/hg38.conservedElts.v3.bed",
                    SD+"/SegDups/tfbs.hg38.bed" ]

annTables = {t.split("/")[-1]: t for t in annotationTables}


if "method" not in config:
    config["method"] = "Integrated"
rule all:
    input:
        svBed="integrated.bed",
        svRenamedBed="integrated_renamed.bed",
        svCalls=expand("integrated.{svtype}.bed",svtype=svTypes),
        pbInclusive="calls.pb.bed",
        pbcalls=expand("pbfinal.{svtype}.bed",svtype=svTypes),
        pbfinalquery=expand("pbfinal_query.{svtype}.bed",svtype=svTypes),
        pbByType=expand("pbcalls.{svtype}.bed", svtype=svTypes),
        svbnIsect=expand("bn_int_query_overlap.{svtype}.bed",svtype=svTypes),
        svQuery=expand("query.{svtype}.bed",svtype=svTypes),
        bnQuery=expand("bnquery.{svtype}.bed",svtype=svTypes),
        bnBed="bncalls.bed",
        bnCalls=expand("bncalls.{svtype}.bed",svtype=svTypes),
        bnPbOvp=expand("bnquery_ovp.{svtype}.bed",svtype=svTypes),
        bnPbBest=expand("bnquery_best.{svtype}.tab",svtype=svTypes),
        bpQuery=expand("bpquery.{svtype}.bed",svtype=svTypes),        
        svBest=expand("pbbest.{svtype}.tab",svtype=svTypes),        
        svIsect=expand("query_overlap.{svtype}.bed",svtype=svTypes),
        bnBest=expand("bnbest.{svtype}.tab",svtype=svTypes),
        readIsect=expand("read_overlap.{svtype}.tab",svtype=svTypes),
        bnReadIsect=expand("bn_read_overlap.{svtype}.bed",svtype=svTypes),
        svCoverage=expand("bn_read_cov.{svtype}.bed",svtype=svTypes),
        intSvCoverage=expand("int_read_cov.{svtype}.bed",svtype=svTypes),
        ilmnMaster=expand("int.{svtype}.bed",svtype=svTypes),
        ilmnRawSup=expand("int.{svtype}.support",svtype=svTypes),
        intPass=expand("int_pass.{svtype}.tab",svtype=svTypes),
        pcel=expand(config["sample"] + ".{svtype}.{score}.pcEL.tsv",svtype=svTypes, score=["20","100"]),
        intFiltPbFinalIsect=expand("int_pb_isect.{svtype}.bed",svtype=svTypes),
        bnIntersect=expand("bn_pb_isect.{svtype}.bed",svtype=svTypes),
        intbn=expand("int_filt_bn_isect.{svtype}.bed",svtype=svTypes),
        pbfinal_bn_annot=expand("pbfinal_bn_annot.{svtype}.tab",svtype=svTypes),
        bn_filt=expand("bn_filt.{svtype}.bed",svtype=svTypes),
        pborthotag=expand("pbfinal_ortho_annot.{svtype}.tab",svtype=svTypes),
        intorthotag=expand("int_filt_ortho_annot.{svtype}.tab",svtype=svTypes),
        bntag=expand("bn_filt_ortho_annot.{svtype}.tab",svtype=svTypes),
        bnorthotag=expand("pbfinal_ortho_annot.{svtype}.tab",svtype=svTypes),
        intfinalpbisect=expand("intfinal_query_pb_isect.{svtype}.bed",svtype=svTypes),
        pbfinal_int_tag=expand("pbfinal_int_tag.{svtype}.tab",svtype=svTypes),
        pbtagbed=expand("pbfinal_ortho.{svtype}.bed",svtype=svTypes),
        bntagbed=expand("bnfinal_ortho.{svtype}.bed",svtype=svTypes),
        inttagbed=expand("intfinal_ortho.{svtype}.bed",svtype=svTypes),
        inttagbedsummary=expand("intfinal_ortho.{svtype}.bed.pbsummary",svtype=svTypes),        
        mergedortho=expand("merged_ortho.{svtype}.bed",svtype=svTypes),
        mergedintegrated=expand("{sample}.merged_nonredundant.{svtype}.bed", sample=config["sample"], svtype=svTypes),
        mnrLocus=expand("{sample}.merged_nonredundant_loci.DEL.bed", sample=config["sample"]),        
        mergedintegratedvcf=expand("{sample}.merged_nonredundant.{svtype}.vcf", sample=config["sample"], svtype=svTypes),
        mnrLocusVCF=expand("{sample}.merged_nonredundant_loci.DEL.vcf", sample=config["sample"]),
        combinedvcf=expand("{sample}.merged_nonredundant.vcf", sample=config["sample"]),        
        mergedintegratedvcfexonstr=expand("{sample}.merged_nonredundant_exon.{svtype}.bed", sample=config["sample"], svtype=svTypes),
        mergedintegratedvcfexonsmerged=expand("{sample}.merged_nonredundant_exon.{svtype}.bed.tr", sample=config["sample"], svtype=svTypes),
        mergedintegratedvcfexons=expand("{sample}.merged_nonredundant_exon.{svtype}.bed.merged", sample=config["sample"], svtype=svTypes),                
        mergedintegratedvcfexonsrvis=expand("{sample}.merged_nonredundant_exon.{svtype}.bed.rvis", sample=config["sample"], svtype=svTypes),
        mergedintegratedvcfexonsrvistab=expand("{sample}.merged_nonredundant_exon.{svtype}.bed.rvis.tab", sample=config["sample"], svtype=svTypes),                
        decorateexons=expand("{sample}.exon_decorate.{svtype}.txt", sample=config["sample"], svtype=svTypes),
        decorateexonsSummary=expand("{sample}.exon_decorate.{svtype}.txt.summary", sample=config["sample"], svtype=svTypes),        
        mergedIntegratedSummary=expand("{sample}.merged_nonredundant.{svtype}.bed.summary", sample=config["sample"], svtype=svTypes),
        mergedIntegratedReclassified=expand("{sample}.merged_nonredundant.{svtype}.bed.rec", sample=config["sample"], svtype=svTypes),        
        exonSummary=expand("{sample}.merged_nonredundant_exon_summary.{svtype}.bed",sample=config["sample"], svtype=svTypes),
        barplots=expand("IntegratedBarchart.{sample}.{svtype}.pdf",sample=config["sample"],svtype=svTypes),
        passpdf=expand("{sample}.{svtype}.pass.pdf",sample=config["sample"], svtype=svTypes),
        passfail=expand("integrated.{svtype}.bed.passfail",svtype=svTypes),
        passmethodpdf=expand("{sample}.{svtype}.method_count.pdf",sample=config["sample"],svtype=svTypes),

        callers=expand("callers.{svtype}.bed", svtype=svTypes),
        callerstab=expand("callers.{svtype}.tab", svtype=svTypes),
        methodpca=expand("MethodPCA.{svtype}.{sample}.pdf",svtype=svTypes, sample=config["sample"]),
        jaccard=expand("Jaccard.{svtype}.{sample}.pdf",svtype=svTypes, sample=config["sample"]),
        methodsbar=expand("MethodsBar.{svtype}.{sample}.pdf",svtype=svTypes, sample=config["sample"]),
        tab=expand("MethodSummary.{svtype}.{sample}.tsv",svtype=svTypes, sample=config["sample"]),
        orthset=expand("int_orthset.{svtype}.bed",svtype=svTypes),
        orthsetquery=expand("int_orthset_query.{svtype}.tab",svtype=svTypes),        
        intsetquery=expand("int_orthset_mapped.{svtype}.bed",svtype=svTypes),
        intsettable=expand("int_orthset_table.{svtype}.bed",svtype=svTypes),    
        annot=expand("int_orthset_table_annotated.{svtype}.bed",svtype=svTypes),
        plotsmall=expand("MethodRecallSmall.{sample}.{svtype}.pdf",sample=config["sample"],svtype=svTypes),
        plotlarge=expand("MethodRecallLarge.{sample}.{svtype}.pdf",sample=config["sample"],svtype=svTypes),
        filt_caller_full=expand("int_caller_full.{svtype}.bed",svtype=svTypes),
        filt_caller=expand("int_caller.{svtype}.bed",svtype=svTypes),
        sensUnion22=expand("SensitivitySpecificity.{sample}.{svtype}.2.2.pdf",sample=config["sample"],svtype=svTypes),
        sensUnion=expand("SensitivitySpecificity.{sample}.{svtype}.2.1.pdf",sample=config["sample"],svtype=svTypes),        
        sensIsect=expand("SensitivitySpecificity.{sample}.{svtype}.3.2.pdf",sample=config["sample"],svtype=svTypes),
        sensTable21=expand("SensitivitySpecificity.{sample}.{svtype}.2.2.tsv",sample=config["sample"],svtype=svTypes),
        sensTable20=expand("SensitivitySpecificity.{sample}.{svtype}.2.1.tsv",sample=config["sample"],svtype=svTypes),
        sensTable32=expand("SensitivitySpecificity.{sample}.{svtype}.3.2.tsv",sample=config["sample"],svtype=svTypes),
        integrated_pbfinal=expand("integrated_pbfinal_isect.{svtype}.bed",svtype=svTypes),
        classplot=expand("{sample}.{svtype}.exon_class.pdf",sample=config["sample"],svtype=svTypes),
        classplottsv=expand("{sample}.{svtype}.exon_class.tsv",sample=config["sample"],svtype=svTypes),
        tabann=expand("annotate-{sample}.{svtype}_table_{table}.tsv",sample=config["sample"],svtype=svTypes, table=annTables.keys()),
        tabannsummary=expand("annotate-{sample}.{svtype}_table_{table}.tsv.union-summary",sample=config["sample"],svtype=svTypes, table=annTables.keys()),
        pbb="pbonly." + config["sample"]+".bed",
        pbreport="pbonly." + config["sample"]+".report.txt",
        pbv="pbonly." + config["sample"]+".vcf",
        plotpbonly="Summary."+config["sample"]+".pbonly.large.pdf",
        tablepbonly="Summary."+config["sample"]+".pbonly.DEL.tsv",
        ilonly=expand("{sample}.il-only.{svtype}.bed",sample=config["sample"], svtype=svTypes),
        ilonlygenes=expand("{sample}.il-only-exons.{svtype}.bed",sample=config["sample"], svtype=svTypes),
        intkv=expand("int_filt_kv.{svtype}.tab", svtype=svTypes),
        intval=expand("int_filt_value.{svtype}.tab", svtype=svTypes),
        intvalbed=expand("int_filt_value.{svtype}.bed", svtype=svTypes),
        reverse=expand("intfinal_reverse.{svtype}.tab",svtype=svTypes),
        scoredSensTable=expand("SensitivitySpecificity.{sample}.{svtype}.{thresh}.scores.tsv",sample=config["sample"], svtype=svTypes,thresh=["2.1", "3.2"]),
        pli=expand("{sample}.merged_nonredundant_exon.{svtype}.bed.pli",sample=config["sample"], svtype=svTypes),
        decorated=expand("{s}.properties_of_unified_callset.{svtype}.bed", s=config["sample"],svtype=svTypes),
        excludeSimpleAndMEI=expand("excluding_simple_and_mei.{op}.bed",op=svTypes),
        netSimplePB=expand("excluding_simple_and_mei_pb_gain.{op}.tab",  op=svTypes),
        UniqueExons=expand("excluding_simple_and_mei_exon.{op}.bed", op=svTypes),
        GainAnnotationSummary=expand("Gain_annotation_Summary.{op}.tab",op=svTypes),
        simpleReciprocal=expand("simple_reciprocal.{op}.bed",op=svTypes),
        browsertrack=expand(config["sample"]+".{b}",b=["bed9","bb"]),
        addAlgorithm=expand("mnr_with_algo.{op}.bed",op=svTypes)
        
        
        
        
        

localBams = { "0": config["localasm"][0], "1": config["localasm"][1] }
if "bn_read_overlap" not in config:
    config["bn_read_overlap"] = "NONE"

if "bn_read_cov" not in config:
    config["bn_read_cov"] = "NONE"


if "softclip" not in config:
    config["softclip"] = "NONE"

rule AddAlg:
    input:
       svs=config["sample"]+".merged_nonredundant.{op}.bed",
       intIsect="intfinal_query_pb_isect.{op}.bed",
       intFile="int_filt.{op}.bed"
    output:
       mnr="mnr_with_algo.{op}.bed"
    params:
       sd=SD,
    shell:"""
{params.sd}/AddAlgorithm.py <(paste {input.intIsect}  <( bioawk -c hdr '{{ print $ALGORITHM;}}' < {input.intFile} ) ) {input.svs} int_filt_pb_key | \
    bioawk -c hdr '{{ print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$svLen"\\t"$svAnn"\\t"$union"\\t"$svAnn"\\t"$svClass"\\t"$ALGORITHM;}}'  > {output.mnr}
"""

rule MakeBB:
    input:
        svs=expand(config["sample"]+".merged_nonredundant.{op}.bed",op=svTypes),
    output:
        bed9=config["sample"]+".bed9",
        bb=config["sample"]+".bb",
    params:
        ref=config["ref"]
    shell:"""
cat {input.svs} | bioawk -c hdr '{{ \
 if ($svType == "deletion") {{ \
   if ($union=="PacBio") {{ rgb="75,0,0";}} \
   else if ($union=="Illumina") {{ rgb="0,75,0";}} \
   else if ($union=="BioNano") {{ rgb="0,0,75";}} \
   else if ($union=="PacBio,Illumina") {{ rgb="75,75,0";}} \
   else if ($union=="PacBio,BioNano") {{ rgb="75,0,75";}} \
   else if ($union=="Illumina,BioNano") {{ rgb="0,75,75";}} \
   else if ($union=="All") {{ rgb="75,75,75";}} }}
 else {{ \
   if ($union=="PacBio") {{ rgb="255,0,0";}} \
   else if ($union=="Illumina") {{ rgb="0,255,0";}} \
   else if ($union=="BioNano") {{ rgb="0,0,255";}} \
   else if ($union=="PacBio,Illumina") {{ rgb="255,255,0";}} \
   else if ($union=="PacBio,BioNano") {{ rgb="255,0,255";}} \
   else if ($union=="Illumina,BioNano") {{ rgb="0,255,255";}} \
   else if ($union=="All") {{ rgb="255,255,255";}} }} \
 print $1"\\t"$2"\\t"$3"\\t"$svType"-"$union"-"$hap"\\t1000\\t+\\t"$2"\\t"$3"\\t"rgb;}}' | bedtools sort > {output.bed9}
 bedToBigBed {output.bed9} {params.ref}.fai {output.bb} -type=bed9
"""

rule SimpleReciprocal:
    input:
        il="integrated.{op}.bed",
        pb="pbfinal.{op}.bed"
    output:
        sr="simple_reciprocal.{op}.bed"
    shell:"""
nf=`head -1 {input.il} | awk '{{ print NF;}}'`
fld=$(($nf+1))
bedtools intersect -header -loj -f 0.5 -r  -a {input.il} -b {input.pb} | \
  bedtools groupby -header -g 1-5 -c $nf -o first -full | \
  bioawk -c hdr '{{ if ($(NF-1) != ".") {{ n+=1;l+=$svLen;}} }} END{{ print NR"\\t"n;}}' > {output.sr}
"""
    
rule CountIlluminaRawSupport:
    input:
        tab="integrated.{svtype}.bed"
    output:
        ilmnRawSup="int.{svtype}.support",
    params:
        fofn=config["fofn"],
        sd=SD,
    shell:"""
{params.sd}/CheckReadAgreementOfSV.py --fofn {params.fofn}  --bed {input.tab} --outFile {output.ilmnRawSup}
"""

rule MakeGainAnnotationSummary:
    input:
        comb=config["sample"] + ".merged_nonredundant.{op}.bed",
    output:
        gain="Gain_annotation_Summary.{op}.tab",
    params:
        sd=SD,
    shell:"""
# 1. Get totals
cat {input.comb} | bioawk -c hdr '{{ if ($union == "PacBio" || $union == "PacBio,BioNano") {{ n+=1; nBP+=$3-$2;}} }} END{{ print "Total\\t"n"\\t"nBP;}}' > {output.gain}

#. 2 Count telomeric SVs (general).
bedtools intersect -header -a {input.comb} -b {params.sd}/../../../regions/hg38.telomeres_500k.bed -u | bioawk -c hdr '{{ if ($union == "PacBio" || $union == "PacBio,BioNano") {{ n+=1; nBP+=$3-$2;}} }} END{{ print "Telomeric\\t"n"\\t"nBP;}}' >> {output.gain}

#. Tandem repeats.
cat {input.comb}| \
   bioawk -c hdr '{{ if ( ($union == "PacBio" || $union == "PacBio,BioNano") && (($is_trf == "TR" && ($svAnn == "NONE" || $svAnn == "0"))  || $svAnn == "TandemRepeat" )) \
        {{ n+=1; nBP+=$3-$2;}} }} END{{ print "Tandem Repeat\\t"n"\\t"nBP;}}'>> {output.gain}

#. Unique 
       
cat {input.comb}| \
   bioawk -c hdr '{{ if ( ($union == "PacBio" || $union == "PacBio,BioNano") && ($is_trf != "TR" && $svAnn != "TandemRepeat" && $svRep < 0.8 )) \
        {{ n+=1; nBP+=$3-$2;}} }} END{{ print "Not-repetitive\\t"n"\\t"nBP;}}'>> {output.gain}

# Mobile elements
cat {input.comb}| \
   bioawk -c hdr '{{ if ( ($union == "PacBio" || $union == "PacBio,BioNano") && \
        ($svAnn != "TandemRepeat") && \
        ($svAnn != "NONE") && \
        (match($svAnn,",") == 0) && \
        ($svAnn != "0") && \
        ($svRep >= 0.8) && \
        (match($svAnn, "Alu") > 0 || match($svAnn, "SVA") > 0 || match($svAnn, "HERV") > 0 || match($svAnn,"L1") > 0)) {{
         n+=1; nBP+=$3-$2;}} }} END{{ print "Mobile element\\t"n"\\t"nBP;}}'>> {output.gain}

# coding
cat {input.comb}| \
   {params.sd}/../../../sv/utils/ToPoint.sh | \
   bedtools intersect -header -a stdin -b {params.sd}/../../../regions/Exons.NoUTR.bed -u | \
   bioawk -c hdr '{{ if ( ($union == "PacBio" || $union == "PacBio,BioNano") && $svAnn !="locus" ) \
        {{ n+=1; nBP+=$svLen;}} }} END{{ print "Coding\\t"n"\\t"nBP;}}' >> {output.gain} 
    
"""

rule MakeSimpleAndNotMEI:
    input:
        comb=config["sample"] + ".merged_nonredundant.{op}.bed"
    output:
        uniq="excluding_simple_and_mei.{op}.bed",
        tab="excluding_simple_and_mei_pb_gain.{op}.tab",
    params:
        sge_opts=config["sge_small"],
        s=config["sample"]
    shell:"""
cat {input.comb}| \
   bioawk -c hdr '{{ if (NR==1 || ($is_trf != "TR" && $svAnn != "TandemRepeat" &&  $svRep< 0.8 )) \
   print $_chrom"\\t"$tStart"\\t"$tEnd"\\t"$hap"\\t"$svType"\\t"$svLen"\\t"$svClass"\\t"$callset"\\t"$union"\\t"$svAnn"\\t"$svRep"\\t"$fracTR"\\t"$is_trf;}}' >{output.uniq}

cat {output.uniq} | bioawk -c hdr '{{ \
  if ($union == "PacBio" || $union  == "PacBio,BioNano") {{ nPB++; lPB+=$svLen; }} \
  if ($union == "PacBio,Illumina") {{ nPBIL++; lPBIL+=$svLen; }} }} \\
  END {{ print {params.s}"\\t"nPB"\\t"lPB"\\t"lPB/nPB"\\t"nPBIL"\\t"lPBIL"\\t"lPBIL/nPBIL;}}' > {output.tab}
        
"""
   
rule MakeSimpleAndNotMEIExons:
    input:
        uniq="excluding_simple_and_mei.{op}.bed",
    output:
        exons="excluding_simple_and_mei_exon.{op}.bed",
    params:
        sge_opts=config["sge_small"],
        sd=SD,
    shell:"""
cat {input.uniq} | {params.sd}/../../../sv/utils/ToPoint.sh  |  bedtools intersect  -a stdin -b {params.sd}/../../../regions/Exons.NoUTR.bed -u -wa| grep -v locus  > {output.exons}
"""

rule ReclassifyMergedNonredundant:
    input:
        mnr=expand("{sample}.merged_nonredundant.{{svtype}}.bed", sample=config["sample"]),
    output:
        rec=expand("{sample}.merged_nonredundant.{{svtype}}.bed.rec", sample=config["sample"]),
    params:
        sd=SNAKEMAKE_DIR,
        sge_opts="-cwd -pe serial 1 -l h_rt=1:00:00 -l mfree=1G"
    shell:"""
{params.sd}/../Enrichment/Reclassify.py --svBed {input.mnr} --out {output.rec}
"""
rule PBOnlyBed:
    input:
        mnr=[config["sample"] + ".merged_nonredundant_loci.DEL.bed", config["sample"] + ".merged_nonredundant.INS.bed"],
    output:
        pbb="pbonly.{sample}.bed"
    params:
        sd=SNAKEMAKE_DIR,
        sge_opts="-cwd -pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        ref=config["ref"]
    shell:"""
cat {input.mnr[0]} | bioawk -c hdr '{{ if ( NR==1 || $union  == "PacBio" || $union == "PacBio,BioNano") print; }}' > {output.pbb}
cat {input.mnr[1]} | bioawk -c hdr '{{ if ( $union  == "PacBio" || $union == "PacBio,BioNano") print; }}' >> {output.pbb}
"""

rule PBOnlyVcf:
    input:
        pbb="pbonly." + config["sample"]+".bed"
    output:
        pbv="pbonly." + config["sample"]+".vcf"   
    params:
        sd=SNAKEMAKE_DIR,
        sge_opts="-cwd -pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        ref=config["ref"],
        sample=config["sample"]
    shell:"""

{SNAKEMAKE_DIR}/../../utils/variants_bed_to_vcf.py --bed {input.pbb} --reference {params.ref} --vcf /dev/stdout  --sample {params.sample} --type sv --addci 5 --fields NALT nAlt NREF nRef SVANN svAnn SVREP svRep SVCLASS svClass NTR nTR BN bnKey SOURCE source --source PHASED-SV --info \"##INFO=<ID=NALT,Number=1,Type=Integer,Description=\\\"Number of reads supporting variant\\\">" "##INFO=<ID=NREF,Number=1,Type=Integer,Description=\\\"Number of reads supporting reference\\\">" "##INFO=<ID=SVANN,Number=1,Type=String,Description=\\\"Repeat annotation of variant\\\">" "##INFO=<ID=SVREP,Number=1,Type=Float,Description=\\\"Fraction of SV annotated as mobile element or tandem repeat\\\">"  "##INFO=<ID=SVCLASS,Number=1,Type=String,Description=\\\"General repeat class of variant\\\">"  "##INFO=<ID=NTR,Number=1,Type=Integer,Description=\\\"Number of tandem repeat bases\\\">"  "##INFO=<ID=BN,Number=1,Type=Float,Description=\\\"Overlap with BioNanoGenomics call.\\\">"  "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\\\"Source method of call, with additional information describing haplotype.\\\">" | bedtools sort -header > {output.pbv}
"""

rule PBOnlyReport:
    input:
        pbv="pbonly." + config["sample"]+".vcf"
    output:
        pbreport="pbonly." + config["sample"]+".report.txt"
    params:
        sd=SNAKEMAKE_DIR,
        sge_opts="-cwd -pe serial 1 -l h_rt=1:00:00 -l mfree=1G"
    shell:"""
{params.sd}/../../PrintMEIReport.py  --vcf {input.pbv} > {output.pbreport}
"""

    
rule IlOnlyGenes:
    input:
        ilonly=config["sample"] + ".il-only.{svtype}.bed",
        exons=SD+"/SegDups/Exons.9.20.bed"
    output:
        ilonlygenes=config["sample"] + ".il-only-exons.{svtype}.bed"
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"],
    shell:"""
cat {input.ilonly} | {params.sd}/Transform.py {wildcards.svtype} | bedtools intersect -a stdin -b {input.exons} -loj | bedtools groupby -g 1-3 -c 7 -o first -full | awk '{{ if ($4 != ".") print; }}' | bedtools groupby -g 4-6 -c 7,10 -o first,first > {output.ilonlygenes}
"""

rule SelectIlOnly:
    input:
        calls=config["sample"] + ".merged_nonredundant.{svtype}.bed"
    output:
        ilonly=config["sample"] + ".il-only.{svtype}.bed"
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"],
    shell:"""
egrep "^#|Illumina" {input.calls} | grep -v "PacBio" | cut -f 1-3 > {output.ilonly}
"""

rule AnnotatePhastConsEL:
    input:
        calls=config["sample"] + ".merged_nonredundant.{svtype}.bed"
    output:
        pcel=config["sample"] + ".{svtype}.{score}.pcEL.tsv"
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"],
    shell:"""
        {params.sd}/../../utils/PrintFractionTR.sh {input.calls} {wildcards.svtype} {params.sd}/SegDups/phastCons{wildcards.score}way.EL.bed --count --header phc_frac phc_count > {output.pcel}
"""

rule PlotPBUnique:
    input:
        calls=expand("{sample}.merged_nonredundant.{svtype}.bed",sample=config["sample"], svtype=svTypes)
    output:
        plot="Summary."+config["sample"]+".pbonly.large.pdf",
        table="Summary."+config["sample"]+".pbonly.DEL.tsv",
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"],
    shell:"""
bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $union == "PacBio" || $union == "Pacbio,BioNano") print ;}}' < {input.calls[0]} > {params.sample}.pbonly.bed
grep -v "^#" {input.calls[1]} >> {params.sample}.pbonly.bed
Rscript {params.sd}/../../utils/PlotStackedBarChart.R --table {params.sample}.pbonly.bed  --reflect T --main "{params.sample} PacBio unique" --outfile Summary.{params.sample}.pbonly --sample {params.sample}

"""

        
rule AnnotatePhastCons:
    input:
        calls=config["sample"] + ".merged_nonredundant.{svtype}.bed",
        pc=lambda wildcards: annTables[wildcards.table],
    output:
        ann="annotate-" + config["sample"] + ".{svtype}_table_{table}.tsv"
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"],
    shell:"""

{params.sd}/../../utils/PrintFractionTR.sh {input.calls} {wildcards.svtype} {input.pc} --count --header phc_frac phc_count > {output.ann}
"""

rule SummarizePhastCons:
    input:
        calls=config["sample"] + ".merged_nonredundant.{svtype}.bed",
        ann="annotate-" + config["sample"] + ".{svtype}_table_{table}.tsv"
    output:
        summary="annotate-" + config["sample"] + ".{svtype}_table_{table}.tsv.union-summary"
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"],
    shell:"""

paste {input.calls} {input.ann} | bioawk -c hdr '{{ if ($phc_count > 0) {{ print $1"\\t"$2"\\t"$3"\\t"$union"\\t"$svLen"\\t"$phc_frac"\\t"$phc_count;}} }} ' | {params.sd}/SummarizeUnion.py > {output.summary}

paste {input.calls} {input.ann} | bioawk -c hdr '{{ if ($phc_count > 0) {{ print $1"\\t"$2"\\t"$3"\\t"$union"\\t"$svLen"\\t"$phc_frac"\\t"$phc_count;}} }} '  > {output.summary}.phc-full
        
"""





      

rule PlotExonCount:
    input:
        exons=config["sample"]+".merged_nonredundant_exon.{svtype}.bed",
    output:
        exonplot=config["sample"]+".{svtype}.exon_class.pdf",
        exontsv=config["sample"]+".{svtype}.exon_class.tsv",
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"],
    shell:"""
Rscript {params.sd}/TabulateExons.R --tab {input.exons} --sample {params.sample} --op {wildcards.svtype}
"""

rule CalcNewFracTR:
    input:
        exons=config["sample"]+".merged_nonredundant_exon.{svtype}.bed",
    output:
        exontr=config["sample"]+".merged_nonredundant_exon.{svtype}.bed.tr",
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"],
    shell:"""
#
# Fraction of sites overlapping tr
#
if [ {wildcards.svtype} == "DEL" ]; then 
cat {input.exons} | {params.sd}/Transform.py {wildcards.svtype} | bedtools intersect -a stdin -b {params.sd}/../../../regions/tandem_repeats_strs.bed -loj | bedtools groupby -g 1-3 -c 16,17 -o collapse,collapse  | {params.sd}/../../utils/PrintFractionTR.py --header frac_tr > {output.exontr}
else
echo "frac_tr" > {output.exontr}
cat {input.exons} | bioawk -c hdr '{{ if (NR > 1) {{b=length($svSeq); print ">"b"#"$svSeq;}} }}' | {params.sd}/../../analysis/RecalcFracTRLine.sh >> {output.exontr}
fi
    
"""

rule DecorateExonsSetdup:
    input:
        exons=config["sample"]+".merged_nonredundant_exon.{svtype}.bed",
        tr=config["sample"]+".merged_nonredundant_exon.{svtype}.bed.tr",        
    output:
        exon_decorate=config["sample"]+".exon_decorate.{svtype}.txt",
        exon_decoratesummary=config["sample"]+".exon_decorate.{svtype}.txt.summary",        
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"],
    shell:"""

#
#  Count fraction of sites overlapping segdup.
#
cat {input.exons} | {params.sd}/Transform.py {wildcards.svtype} | bedtools intersect -a stdin -b {params.sd}/SegDups/grch38_superdups.merged.bed -loj | bedtools groupby -g 1-3 -c 16,17 -o collapse,collapse  | {params.sd}/../../utils/PrintFractionTR.py --header frac_segdup > {params.sample}.exons.{wildcards.svtype}.segdup


paste {params.sample}.exons.{wildcards.svtype}.segdup {input.tr}  > {output.exon_decorate}



#
# Summarize all events.
#

# Segdup
paste {input.exons} {output.exon_decorate} | bioawk -c hdr '{{ if (substr($0,0,1) == "#") {{ print "#union"; }} else {{ if ($frac_segdup > 0.5) {{ print $union;}} }} }}'  | {params.sd}/UnionTable.py --noHeader  > {output.exon_decorate}.summary.1

    # Tandem repeat
paste {input.exons} {output.exon_decorate} | bioawk -c hdr '{{ if (substr($0,0,1) == "#") {{ print "#union"; }} else {{ if ($frac_tr > 0.5) {{ print $union;}} }} }}'  | {params.sd}/UnionTable.py  --noHeader > {output.exon_decorate}.summary.2

# exon segdup
paste {input.exons} {output.exon_decorate} |  egrep "^#|EXON"  |  bioawk -c hdr '{{ if (substr($0,0,1) == "#") {{ print "#union"; }} else {{ if ($frac_segdup > 0.5) {{ print $union;}} }} }}'  | {params.sd}/UnionTable.py --noHeader  > {output.exon_decorate}.summary.3

paste {input.exons} {output.exon_decorate} | bioawk -c hdr '{{ if (substr($0,0,1) == "#") {{ print "#union"; }} else {{ if ($frac_segdup > 0.5) {{ print $union;}} }} }}'  | {params.sd}/UnionTable.py --noHeader  > {output.exon_decorate}.summary.1

paste {input.exons} {output.exon_decorate} | bioawk -c hdr '{{ if (substr($0,0,1) == "#") {{ print "#union"; }} else {{ if ($frac_tr > 0.5) {{ print $union;}} }} }}'  | {params.sd}/UnionTable.py  --noHeader > {output.exon_decorate}.summary.2


paste {input.exons} {output.exon_decorate} |  egrep "^#|EXON"  |  bioawk -c hdr '{{ if (substr($0,0,1) == "#") {{ print "#union"; }} else {{ if ($frac_segdup > 0.5) {{ print $union;}} }} }}'  | {params.sd}/UnionTable.py --noHeader  > {output.exon_decorate}.summary.3

paste {input.exons} {output.exon_decorate} | egrep "^#|EXON" |  bioawk -c hdr '{{ if (substr($0,0,1) == "#") {{ print "#union"; }} else {{ if ($frac_tr > 0.5) {{ print $union;}} }} }}'  | {params.sd}/UnionTable.py  --noHeader > {output.exon_decorate}.summary.4

paste  {output.exon_decorate}.summary.1 {output.exon_decorate}.summary.2 {output.exon_decorate}.summary.3 {output.exon_decorate}.summary.4 > {output.exon_decorate}.summary

rm -f {output.exon_decorate}.summary.1 {output.exon_decorate}.summary.2 {output.exon_decorate}.summary.3 {output.exon_decorate}.summary.4  
    
"""


        
rule MakeExonSummary:
    input:
        exons=config["sample"]+".merged_nonredundant_exon.{svtype}.bed",
    output:
        exonSummary=config["sample"]+".merged_nonredundant_exon_summary.{svtype}.bed",
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
    shell:"""
egrep "^#|EXON" {input.exons} | bedtools groupby -header -g 10 -c 5  -o sum -full | tr -d "()" |  {params.sd}/UnionTable.py > {output.exonSummary}.1
egrep "^#|EXON" {input.exons} |  {params.sd}/UnionTable.py --noHeader > {output.exonSummary}.2
egrep "^#|UTR" {input.exons} | bedtools groupby -header -g 10 -c 5  -o sum -full | tr -d "()" | {params.sd}/UnionTable.py  --noHeader > {output.exonSummary}.3
egrep "^#|UTR" {input.exons} | {params.sd}/UnionTable.py  --noHeader > {output.exonSummary}.4

paste {output.exonSummary}.1 {output.exonSummary}.2 {output.exonSummary}.3 {output.exonSummary}.4 > {output.exonSummary}

/bin/rm -f paste {output.exonSummary}.1 {output.exonSummary}.2 {output.exonSummary}.3 {output.exonSummary}.4         
"""

rule MakeDecoratedExonSummary:
    input:
        exonSummary=config["sample"]+".merged_nonredundant_exon_summary.{svtype}.bed",
        exon_decorate=config["sample"]+".exon_decorate.{svtype}.txt.summary",
    output:
        decorated=config["sample"]+".properties_of_unified_callset.{svtype}.bed",
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
    shell:"""
paste {input.exonSummary} {input.exon_decorate} > {output.decorated}
"""
        
rule MakeIntOrthTable:
    input:
        query="int_orthset_mapped.{svtype}.bed"
    output:
        table="int_orthset_table.{svtype}.bed"
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
    shell:"""
{params.sd}/CallersToMarginTable.py  --svbed {input.query} --callerKey ALGORITHM_db > {output.table}
"""

rule MakeIntCallerTable:
    input:
        filt="int.{svtype}.bed"
    output:
        filt_caller="int_caller.{svtype}.bed"
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
    shell:"""
{params.sd}/CallersToMarginTable.py  --svbed {input.filt} --out {output.filt_caller} --callerKey ALGORITHM
"""

rule MakeIntFullCallerTable:
    input:
        filt="int.{svtype}.bed",
        filt_caller="int_caller.{svtype}.bed"
    output:
        filt_caller_full="int_caller_full.{svtype}.bed"
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
    shell:"""
paste {input.filt} {input.filt_caller} > {output.filt_caller_full}
"""

    
rule MakeIntOrthAnnotatedTable:
    input:
        table="int_orthset_table.{svtype}.bed",
        orthset="int_orthset.{svtype}.bed"
    output:
        annot="int_orthset_table_annotated.{svtype}.bed"
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
    shell:"""
bioawk -c hdr '{{ print $_chrom"\\t"$tStart"\\t"$tEnd"\\t"$svLen"\\t"$union;}}' < {input.orthset} > {output.annot}.tmp
paste {output.annot}.tmp {input.table} > {output.annot}
"""

rule PlotIntOrthAnnotatedTable:
    input:
        annot="int_orthset_table_annotated.{svtype}.bed",
        filt="int_caller_full.{svtype}.bed"
    output:
        plotsmall="MethodRecallSmall.{sample}.{svtype}.pdf",
        plotlarge="MethodRecallLarge.{sample}.{svtype}.pdf"
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"]
    shell:"""
Rscript {params.sd}/../../../plotting/plot_method_recall.R --tab {input.annot} --filt {input.filt} --sample {params.sample} --op {wildcards.svtype}
"""
    
rule MakeIntOrthSet:
    input:
        nonredundant=config["sample"]+".merged_nonredundant.{svtype}.bed"
    output:
        orthset="int_orthset.{svtype}.bed"
    params:
        sge_opts=config["sge_small"],    
    shell:"""
bioawk -c hdr '{{ if (substr($0,0,1) == "#") {{ print;}} else {{ if ($union == "PacBio,Illumina" || $union == "PacBio" || $union == "BioNano"|| $union == "Pacbio,BioNano" || $union == "BioNano,Illumina" || $union == "ALL") print;}} }}' <  {input.nonredundant} > {output.orthset}
"""

rule MakeIntQuery:
    input:
        integrated="integrated.{svtype}.bed",
        orthset="int_orthset.{svtype}.bed"
    output:
        orthsetquery="int_orthset_query.{svtype}.tab"
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR,
    shell:"""
cat {input.orthset} | \
awk '{{ print $1"\\t"$2"\\t"$3"\\t"$0;}}' | \
bedtools slop -header  -g {params.ref}.fai -b 1000 -i stdin | \
bedtools intersect -a stdin -b {input.integrated} -header -loj | \
bioawk -c hdr '{{ if (substr($0,0,1) == "#") \
 {{ print $0"\\t#oChrom\\toStart\\toEnd\\tosvType\\tosvLensv\\toSeq\\tCIEND\\tCIPOS\\tEND\\toSV_TYPE\\toALGORITHM";}} \
  else {{ print;}} }}' | \
 {params.sd}/SelectBestcall.py --header "integrated" --key --qs 4 --qe 5  > {output.orthsetquery}
"""

rule MakeIntToFullTable:
    input:
       orthsetquery="int_orthset_query.{svtype}.tab",
       integrated="integrated.{svtype}.bed"
    output:
       intsetquery="int_orthset_mapped.{svtype}.bed"
    params:       
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
{params.sd}/SelectRowByKey.py --keyFile {input.orthsetquery} --keyColumn integrated_key --dbFile {input.integrated}  > {output.intsetquery}
"""
    
rule MakeCallers:
    input:
        vcf=config["svvcf"]
    output:
        bed="callers.{svtype}.bed"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
{params.sd}/../../utils/variants_vcf_to_bed.py --vcf {input.vcf} --groupCommas --fields ALGORITHM SVLEN CIPOS CIEND  --filter --out /dev/stdout | \
 sed "s/ALU/INS/" | sed "s/ALU;INS/INS/" | sed "s/DUP/INS/" | sed "s/LINE1/INS/" | sed "s/SVA/INS/" | sed "s/TANDUP/INS/" | sed "s/TAN/INS/" | \
egrep "^#|{wildcards.svtype}" > {output.bed}
"""


rule MakeSVTab:
    input:
        callerbed="callers.{svtype}.bed"
    output:
        callertab="callers.{svtype}.tab"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
{params.sd}/CallersToMarginTable.py --callerKey ALGORITHM --svbed {input.callerbed} --out {output.callertab}
"""

rule PlotMethods:
    input:
        callertab="callers.{svtype}.tab",
        filt="int_caller_full.{svtype}.bed"
    output:
        pca="MethodPCA.{svtype}.{sample}.pdf",
        jac="Jaccard.{svtype}.{sample}.pdf",
        meth="MethodsBar.{svtype}.{sample}.pdf",
        tab="MethodSummary.{svtype}.{sample}.tsv"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"]
    shell:"""
Rscript {params.sd}/individual_caller.R --sv {input.filt} --count {input.callertab}  --operation {wildcards.svtype} --sample {params.sample}
"""

ncrlim = {"DEL": "0.10", "INS" : "0.20"}

rule ScoreRuns:
    input:
        sensTable=expand("SensitivitySpecificity.{sample}.{{svtype}}.{{count}}.{{minn}}.tsv",sample=config["sample"]),
        intfilt="int_filt.{svtype}.bed",
        intfinalpbisect="intfinal_query_pb_isect.{svtype}.bed",        
        mnr=expand("{sample}.merged_nonredundant.{{svtype}}.bed.rec",sample=config["sample"])
    output:
        scoredSensTable=expand("SensitivitySpecificity.{sample}.{{svtype}}.{{count}}.{{minn}}.scores.tsv",sample=config["sample"])
    params:
        ncr=lambda wildcards : ncrlim[wildcards.svtype],
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"]
    shell:"""

paste {input.intfilt} {input.intfinalpbisect} > {input.intfilt}.final

rm -f SensitivitySpecificity.{params.sample}.{wildcards.svtype}.{wildcards.count}.{wildcards.minn}.scores.tsv
touch SensitivitySpecificity.{params.sample}.{wildcards.svtype}.{wildcards.count}.{wildcards.minn}.scores.tsv
for comb in `bioawk -c hdr -v ncrlim={params.ncr} '{{ if ($ncr < ncrlim) print; }}' < {input.sensTable} | sort -k1,1nr | grep -v liWGS |  cut -f 4`; do

echo "{params.sd}/ScoreSubset.sh --comb $comb --mincount {wildcards.minn} --sample {params.sample} --operation {wildcards.svtype} --integrated {input.intfilt}.final >> SensitivitySpecificity.{params.sample}.{wildcards.svtype}.{wildcards.count}.{wildcards.minn}.scores.tsv" > /dev/stderr
        
{params.sd}/ScoreSubset.sh --comb $comb --mincount {wildcards.minn} --sample {params.sample} --operation {wildcards.svtype} --integrated {input.intfilt}.final >> SensitivitySpecificity.{params.sample}.{wildcards.svtype}.{wildcards.count}.{wildcards.minn}.scores.tsv
        

  
 done
"""    
    
rule PlotSensSpec:
    input:
        filt="int_caller_full.{svtype}.bed"
    output:
        sensUnion="SensitivitySpecificity.{sample}.{svtype}.3.2.pdf",
        sensIsect="SensitivitySpecificity.{sample}.{svtype}.2.2.pdf",
        sensUnionTwo="SensitivitySpecificity.{sample}.{svtype}.2.1.pdf",    
        sensTable="SensitivitySpecificity.{sample}.{svtype}.2.2.tsv",
        sensTable21="SensitivitySpecificity.{sample}.{svtype}.2.1.tsv",        
        sensTable32="SensitivitySpecificity.{sample}.{svtype}.3.2.tsv",        
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"]
    shell:"""
Rscript {params.sd}/../../../plotting/plot_method_sensitivity.R --tab {input.filt}  --operation {wildcards.svtype} --sample {params.sample}
"""

rule PlotSVs:
    input:
        merged="merged_ortho.{op}.bed"
    output:
        plot="IntegratedBarchart."+config["sample"]+".{op}.pdf"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"],
        method=config["method"]
    shell:"""
Rscript {params.sd}/RenderStackedBarchart.R --tab {input.merged} --sample {params.sample} --op {wildcards.op} --method {params.method}
"""

rule SummarizeOrthoIntegrated:
    input:
        mergedint=config["sample"]+".merged_nonredundant.{op}.bed"
    output:
        summary=config["sample"]+".merged_nonredundant.{op}.bed.summary"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
        sample=config["sample"]
    shell:"""
Rscript {params.sd}/SummarizeMergedNonRedundant.R --tab {input.mergedint} --operation {wildcards.op} --sample {params.sample}
"""

#
#  merged_ortho was formed by first filtering the IL and BN tables, then
# adding dataset membership information, then merging to a fixed number of columns.
#

rule MergeOrthoToIntegrated:
    input:
        merged="merged_ortho.{op}.bed"
    output:
        mergedint=config["sample"]+".merged_nonredundant.{op}.bed"
    params:
        sge_opts=config["sge_small"],
    shell:"""

cat {input.merged} | bioawk -c hdr '{{ if (substr($0,0,1) == "#") {{ print $0"\\tunion";}}\
  else {{ \
     union="NONE"; \
     doPrint="FALSE"; \
     if ($in_pb == "TRUE" && $in_bn == "FALSE" && $in_int == "FALSE" && $callset == "pacbio") {{ union="PacBio"; doPrint="TRUE";}} \
     if ($in_pb == "FALSE" && $in_bn == "TRUE" && $in_int == "FALSE" && $callset == "bionano") {{ union="BioNano"; doPrint="TRUE";}}\
     if ($in_pb == "FALSE" && $in_bn == "FALSE" && $in_int == "TRUE" && $callset == "Illumina") {{ union="Illumina"; doPrint="TRUE";}} \
     if ($in_pb == "TRUE" && $in_bn == "TRUE"  && $in_int == "FALSE" && $callset == "pacbio") {{ union="PacBio,BioNano"; doPrint="TRUE";}} \
     if ($in_pb == "TRUE" && $in_bn == "FALSE"  && $in_int == "TRUE" && $callset == "pacbio") {{ union="PacBio,Illumina"; doPrint="TRUE";}} \
     if ($in_pb == "FALSE" && $in_bn == "TRUE" && $in_int == "TRUE" && $callset == "Illumina") {{ union="Illumina,BioNano"; doPrint="TRUE";}} \
     if ($in_pb == "TRUE" && $in_bn == "TRUE" && $in_int == "TRUE" && $callset == "pacbio") {{ union="All"; doPrint="TRUE";}} \
     if (doPrint == "TRUE") {{ print $0"\\t"union; }} \
   }} }} ' > {output.mergedint}
"""

rule AddLoci:
    input:
        mnr=config["sample"]+".merged_nonredundant.DEL.bed",
        orig=config["pbfinal"]
    output:
        mnrLocus=config["sample"]+".merged_nonredundant_loci.DEL.bed",
    params:
        sge_opts=config["sge_small"],
        sd=SD,
    shell:"""
cat {input.orig} | bioawk -c hdr '{{ if (NR == 1 || $svType == "locus") print;}}' > {output.mnrLocus}.tmp.locus
{params.sd}/../../utils/MergeFiles.py --files {input.mnr} {output.mnrLocus}.tmp.locus --out {output.mnrLocus}
rm -f {output.mnrLocus}.tmp.locus
"""
    
#
# This takes all of the SV tables that have been labeled whether or not
# a call exists in another dataset, and merges them into one master table that
# has a fixed number of columns.  This master table can then be condensed
# to nonredundant calls.
#

rule MergeOrtho:
    input:
        ortho=expand("{meth}final_ortho.{{op}}.bed", meth=["pb","bn","int"])
    output:
        merged="merged_ortho.{op}.bed"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
{params.sd}/../../utils/MergeFiles.py --files {input.ortho} --out {output.merged}
"""
    

rule MakeBNCalls:
    input:
        bnvcf=config["bionano-vcf"]
    output:
        bncalls="bncalls.bed"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
zcat {input.bnvcf} | {params.sd}/../../utils/variants_vcf_to_bed.py --vcf /dev/stdin --out {output.bncalls} --bionano 
"""

rule FilterBNCalls:
    input:
        bnCalls="bncalls.{op}.bed",
        pbQuery="bnquery_best.{op}.tab",
        readOvp="bn_read_overlap.{op}.bed",
        readCov="bn_read_cov.{op}.bed",
        intOvp ="int_filt_bn_isect.{op}.bed"
    output:
        bnFile="bn_filt.{op}.bed"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
    shell:"""
paste {input.bnCalls} {input.pbQuery} {input.readOvp} {input.readCov} {input.intOvp} | \
  bioawk -c hdr -v op="{wildcards.op}" '{{ if (NR==1 || \
          ($pbbest > 0.25 || ($overlap > 0.5 && $number > 5) || op=="INS" || (op=="DEL" && $coverage < 30))) print; }}' > {output.bnFile}
"""
    
        
rule SeparateBNCalls:
    input:
        bncalls="bncalls.bed"
    output:
        bnCallsOp="bncalls.{op}.bed"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,
        operation=lambda wildcards: pbSVTypeMap[wildcards.op]
    shell:"""
nf=`head -1 {input.bncalls} | awk '{{ print NF;}}'`
egrep "^#|{params.operation}" {input.bncalls} | \
  awk '{{ if (substr($0,0,1) == "#" || $3-$2 < 6000000) print; }}' | \
  bedtools groupby -header -g 1-3 -c 4 -o first -full | cut -f 1-$nf | \
  {params.sd}/../../utils/rmdup.py > {output.bnCallsOp}
"""

#
# Merge the support from the pacbio calls, bionano overlap, and read overlap into one.
#
rule IntegratedFilter:
    input:
        intPbSup="pbbest.{svtype}.tab",
        intBnSup="bnbest.{svtype}.tab",
        readOvp="read_overlap.{svtype}.tab",
        intcov="int_read_cov.{svtype}.bed",
        sup="int.{svtype}.support"        
    output:
        intPass="int_pass.{svtype}.tab",        
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        svtype=lambda wildcards: wildcards.svtype
    shell:"""
if [ {params.svtype} == "DEL" ]; then 
   paste {input.intPbSup} {input.intBnSup} {input.readOvp} {input.intcov} {input.sup} | bioawk -c hdr '{{ if (NR==1) {{ print "orth_filter";}} else {{ if (($number >= 3) || $pbbest > 0.25 || $bnoverlap > 0.25 || $coverage < 25 || $nSupport > 3) {{print "PASS";}} else {{ print "FAIL";}} }} }}' > {output.intPass}
else
   paste {input.intPbSup} {input.intBnSup} {input.readOvp} {input.intcov} {input.sup} | bioawk -c hdr '{{ if (NR==1) {{ print "orth_filter";}} else {{ if (($number >= 3) || $pbbest > 0.25 || $bnoverlap > 0.25 || $nSupport > 3) {{print "PASS";}} else {{ print "FAIL";}} }} }}' > {output.intPass}
fi
"""

rule AnnotateIntegratedPass:
    input:
        passtab="int_pass.{svtype}.tab",
        svbed="integrated.{svtype}.bed"
    output:
        passfail="integrated.{svtype}.bed.passfail"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,
        sample=config["sample"]
    shell:"""
paste {input.svbed} {input.passtab} > {input.svbed}.passfail        
"""

rule PlotIlluminaConcordant:
    input:
        passfail="integrated.{svtype}.bed.passfail"    
    output:
        passpdf=config["sample"] + ".{svtype}.method_count.pdf"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,
        sample=config["sample"]
    shell:"""
Rscript {params.sd}/PlotIlluminaPassByMethod.R  --tab {input.passfail}  --operation {wildcards.svtype} --sample {params.sample}
"""

rule PlotIlluminaResolved:
    input:
        passfail="integrated.{svtype}.bed.passfail"    
    output:
        passpdf=config["sample"] + ".{svtype}.pass.pdf"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,
        sample=config["sample"]
    shell:"""
Rscript {params.sd}/PlotIlluminaPass.R --tab {input.passfail}  --operation {wildcards.svtype} --sample {params.sample}
"""

#
# Intesect Illumina calls with pacbio, with an extra annotation that
# states if the intersection is anywhere near the pacbio breakpoints
# (PB_OVP),  not near but overlapping (PB_REFINE), or not touching
# UNIQUE
#

rule IntersectIlluminaFiltWithPbFinal:
    input:
        intFilt="int_filt.{svtype}.bed",
        pbFinal="pbfinal.{svtype}.bed"
    output:
        isect="int_pb_isect.{svtype}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,        
    shell:"""
  echo -e "intdist\\tkey" > {output.isect}
        {{ {{ head -1 {input.intFilt}; head -1 {input.pbFinal} | {params.sd}/../../utils/AddIndex.py _2; echo "\tintdist";}} | tr "\\n" "\\t" ; \
    {{ bedtools sort -i {input.intFilt} | bedtools closest -a stdin -b {input.pbFinal} -d -t first ; }} ; }}  | \
    bioawk -c hdr '{{ minSVLen=$svLen; maxSVLen=$svLen_2; \
                    ovp="NONE"; \
                    key="NONE"; \
                    if ($svLen_2 < $svLen) {{ minSVLen=$svLen_2; maxSVLen=$svLen;}} \
                    if ($intdist < 1000) {{ \
                       key=$_chrom_2"_"$tStart_2"_"$tEnd_2; \
                       if (minSVLen/maxSVLen > 0.5) {{ ovp="PB_OVP";  }} \
                       else {{ ovp="PB_REFINE"; }} \
                    }} \
                    else {{ ovp="UNIQUE"; }} \
                  print ovp"\\t"key; \
                  }}'  >> {output.isect}
"""


rule IntersectBioNanoWithPbFinal:
    input:
        bnfinal="bnquery.{svtype}.bed",
        pbfinal="pbfinal.{svtype}.bed"
    output:
        bnIntersect="bn_pb_isect.{svtype}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
    shell:"""
bedtools intersect -header -a {input.bnfinal} -b {input.pbfinal} -loj | \
 awk '{{ if (substr($0,0,1) == "#") {{ \
   print $0"\\toChrom\\toStart\\toEnd";}} \
        else {{ print $0; }} }}' |\
    ~/projects/HGSVG/hgsvg/sv/analysis/IlluminaCallers/SelectBestcall.py --infile /dev/stdin --ts 14 --te 15 --bnidx  --key --header pbfinal --out {output.bnIntersect}
"""


rule AddBNIsectLabelToPB:
    input:
        bnisect="bn_pb_isect.{svtype}.bed",
        pbfinal="pbfinal.{svtype}.bed"
    output:
        pbfinal_bn_annot="pbfinal_bn_annot.{svtype}.tab"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,
    shell:"""
bioawk -c hdr '{{ print $pbfinal_key;}}' < {input.bnisect}  | {params.sd}/KeyToColumn.py --keys /dev/stdin --bed {input.pbfinal} --header bn_match --out {output.pbfinal_bn_annot}
"""

    
rule MakeIntIsectKV:
    input:
        intrev="intfinal_reverse.{svtype}.tab",
        intfilt="int_filt.{svtype}.bed"
    output:
        intkv="int_filt_kv.{svtype}.tab"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,
    shell:"""
paste {input.intrev} {input.intfilt} | bioawk -c hdr '{{ if (NR ==1) {{ print "key\\tvalue"}} else {{ print $pbkey"\t"$2"_"$3"_"$4; }} }}'  > {output.intkv}
"""


rule MakeIntIsectValue:
    input:
        bed=expand("{sample}.merged_nonredundant.{{svtype}}.bed",sample=config["sample"]),
        intfilt="int_filt_kv.{svtype}.tab"
    output:
        vals="int_filt_value.{svtype}.tab",
        valsbed="int_filt_value.{svtype}.bed"        
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,
    shell:"""
{params.sd}/MakeIntFiltkey.py --mnr {input.bed} --values {input.intfilt} --out {output.vals}
paste {input.bed} {output.vals} > {output.valsbed}
"""



    
rule AddPBCaller:
    input:
        pbfinal="pbfinal.{svtype}.bed"
    output:
        pbtag="pbfinal_tag.{svtype}.tab"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
    shell:"""
echo "in_pacbio" > {output.pbtag}
cat {input.pbfinal} | awk '{{ print "TRUE"; }}' >> {output.pbtag}
"""
    
rule AddBNCaller:
    input:
        bnfinal="bn_filt.{svtype}.bed"
    output:
        bntag="bn_filt_tag.{svtype}.tab"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
    shell:"""
echo "in_bn" > {output.bntag}
cat {input.bnfinal} | awk '{{ print "TRUE"; }}' >> {output.bntag}
"""

rule AddIntCaller:
    input:
        intfinal="int_filt.{svtype}.bed"
    output:
        inttag="int_filt_tag.{svtype}.tab"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
    shell:"""
echo "in_int" > {output.inttag}
cat {input.intfinal} | awk '{{ print "TRUE"; }}' >> {output.inttag}
"""

rule AddPacBioOrthoTags:
    input:
        pbbn="pbfinal_bn_annot.{svtype}.tab",
        pbint="pbfinal_int_tag.{svtype}.tab"
    output:
        tag="pbfinal_ortho_annot.{svtype}.tab"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
    shell:"""
echo -e "callset\tin_pb\tin_bn\tin_int" > {output.tag}
paste {input.pbbn} {input.pbint} | bioawk -c hdr '{{ in_bn="FALSE"; in_int="FALSE";\
   if ($bn_match != "NO") {{ in_bn="TRUE";}} \
   if ($int_match != "NO") {{ in_int="TRUE";}} \
   if ($1 != "bn_match") {{ print "pacbio\tTRUE\t"in_bn"\t"in_int;}} }}' >> {output.tag}
"""

rule AddBNOrthoTags:
    input:
       bnfilt="bn_filt.{svtype}.bed"
    output:
       tag="bn_filt_ortho_annot.{svtype}.tab"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
    shell:"""
echo -e "callset\tin_bn\tin_pb\tin_int" > {output.tag}
cat {input.bnfilt} | bioawk -c hdr '{{ in_pb="FALSE"; in_int="FALSE"; \
  if ($pbbest >= 0.5) {{ in_pb="TRUE"; }} \
  if ($int_isect_key != "NONE") {{ in_int="TRUE"; }} \
  if (substr($0,0,1) != "#") {{ print "bionano\tTRUE\t"in_pb"\t"in_int; }} }}' >> {output.tag}
"""


#
# For the IL callset, record whether or not it is seen in the pacbio
# callsets, given from the external key int_filt_pb_key.
#
rule AddIntOrthoTags:
    input:
       intfilt="int_filt.{svtype}.bed",
       intpb="intfinal_query_pb_isect.{svtype}.bed"
    output:
       tag="int_filt_ortho_annot.{svtype}.tab"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
    shell:"""
echo -e "callset\\tin_int\\tin_pb\\tin_bn" > {output.tag}
paste {input.intfilt} {input.intpb} | bioawk -c hdr '{{ in_pb="FALSE"; in_bn="FALSE"; \
  if ($int_filt_pb_key != "NONE") {{ in_pb="TRUE"; }} \
  if ($bnoverlap >= 0.5) {{ in_bn="TRUE"; }} \
  if (substr($0,0,1) != "#") {{ print "Illumina\tTRUE\t"in_pb"\t"in_bn; }} }}' >> {output.tag}
"""

rule CombinePBOrtho:
    input:
        pbfinal="pbfinal.{svtype}.bed",
        tag="pbfinal_ortho_annot.{svtype}.tab"
    output:
        pbtag="pbfinal_ortho.{svtype}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
    shell:"""
paste {input.pbfinal} {input.tag} > {output.pbtag}
"""

rule CombineBNOrtho:
    input:
        bnFinal="bn_filt.{svtype}.bed",
        tag="bn_filt_ortho_annot.{svtype}.tab"
    output:
        bntag="bnfinal_ortho.{svtype}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
    shell:"""
paste {input.bnFinal} {input.tag} > {output.bntag}
"""

rule CombineIntOrtho:
    input:
        intFinal="int_filt.{svtype}.bed",
        tag="int_filt_ortho_annot.{svtype}.tab"
    output:
        inttag="intfinal_ortho.{svtype}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
    shell:"""
paste {input.intFinal} {input.tag} > {output.inttag}
"""

rule SummarizeIntOrthoPBSupport:
    input:
        inttag="intfinal_ortho.{svtype}.bed"
    output:
        inttagsummary="intfinal_ortho.{svtype}.bed.pbsummary"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
    shell:"""
cat {input.inttag} | bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $in_pb == "FALSE") print;}}' | wc -l > {output.inttagsummary}
cat {input.inttag} | bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $in_pb == "FALSE") print;}}' | bioawk -c hdr '{{ if ($pbbest > 0.5 ) print;}}' | wc -l >> {output.inttagsummary}
"""



rule MakeIntegratedMaster:
    input:
        intOrig="integrated.{svtype}.bed",
        intPbSup="pbbest.{svtype}.tab",
        intBnSup="bnbest.{svtype}.tab",
        readOvp="read_overlap.{svtype}.tab",
        intPass="int_pass.{svtype}.tab",
        intcov="int_read_cov.{svtype}.bed",
        sup="int.{svtype}.support"
    output:
        intMaster="int.{svtype}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
    shell:"""

paste {input.intOrig} {input.intPbSup} {input.intBnSup} {input.readOvp} {input.intPass} {input.intcov} {input.sup} > {output.intMaster}
"""

rule FilteredIntegratedCalls:
    input:
        intMaster="int.{svtype}.bed"
    output:
        intFilt="int_filt.{svtype}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
    shell:"""
bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $orth_filter == "PASS") print;}}' < {input.intMaster} > {output.intFilt}
"""
        
    

rule MakePBFinal:
    input:
        pbfinal=config["pbfinal"]
    output:
        pbfinaltype="pbfinal.{svtype}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,
        pbop=lambda wildcards : pbSVTypeMap[wildcards.svtype]
    shell:"""
egrep "^#|{params.pbop}" {input.pbfinal} > {output.pbfinaltype}
"""

rule MakePBFinalQuery:
    input:
        pbfinal="pbfinal.{svtype}.bed"
    output:
        pbfinalquery="pbfinal_query.{svtype}.bed",
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        ref=config["ref"],
        window=config["window"]
    shell:"""
cat {input.pbfinal} | \
awk '{{ if (substr($0,0,1) == "#") {{ print "#sChrom\\tsStart\\tsEnd\\t"$0;}} else {{ print $1"\\t"$2"\\t"$3"\\t"$0;}} }}'| \
 bedtools  slop -header -i stdin -b {params.window} -g {params.ref}.fai >  {output.pbfinalquery}
"""

#rule IntersectPBFinalQueryAndIntFinal:
#    input:
#        pbfinalquery="pbfinal_query.{svtype}.bed",
#        intfilt="int_filt.{svtype}.bed"
#    output:
#        pbfinalintisect="pbfinal_query_int_isect.{svtype}.bed"
#    params:
#        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
#    shell:"""        
#bedtools intersect -header -a {input.pbfinalquery} -b {input.intfilt} -header -loj | awk '{{ if (substr($0,0,1) == "#") {{ print $0"\\toChrom\\toStart\\toEnd";}} else {{ print $0;}} }}' > {output.pbfinalintisect}
#"""
#


#
# Run sloppy intersection int-> pb to find what elements in the pb
# callset are in the pacbio callset.
#
rule IntersectIntFinalQueryAndPBFinal:
    input:
        pbfinalquery="pbfinal.{svtype}.bed",
        intfilt="int_filt.{svtype}.bed"
    output:
        intfinalpbisect="intfinal_query_pb_isect.{svtype}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        ref=config["ref"],
        window=config["window"],
        sd=SNAKEMAKE_DIR,
    shell:"""
cat {input.intfilt} | \
  awk '{{ if (substr($0,0,1) == "#") {{ print "#sChrom\\tsStart\\tsEnd\\t"$0;}} else {{ print $1"\\t"$2"\\t"$3"\\t"$0;}} }}' | \
  bedtools slop -header -i stdin -b {params.window} -g {params.ref}.fai |\
  bedtools intersect -a stdin -b {input.pbfinalquery} -header -loj | \
  awk '{{ if (substr($0,0,1) == "#") {{ print $0"\\toChrom\\toStart\\toEnd";}} else {{ print $0; }} }} ' | \
  {params.sd}/SelectBestcall.py --infile /dev/stdin --header int_filt_pb  --ts 21 --te 22 --svlen 7  --key > {output.intfinalpbisect}.tmp
paste {output.intfinalpbisect}.tmp <( bioawk -c hdr '{{ print $ALGORITHM; }} < {input.intfilt} ) > {output.intfinalpbisect}
"""

rule IntersectIntegratedOrigAndPBFinal:
    input:
        pbfinalquery="pbfinal.{svtype}.bed",
        integrated="integrated.{svtype}.bed"
    output:
        integrated_pbfinal="integrated_pbfinal_isect.{svtype}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        ref=config["ref"],
        window=config["window"],
        sd=SNAKEMAKE_DIR,
    shell:"""
cat {input.integrated} | \
  awk '{{ if (substr($0,0,1) == "#") {{ print "#sChrom\\tsStart\\tsEnd\\t"$0;}} else {{ print $1"\\t"$2"\\t"$3"\\t"$0;}} }}' | \
  bedtools slop -header -i stdin -b {params.window} -g {params.ref}.fai |\
  bedtools intersect -a stdin -b {input.pbfinalquery} -header -loj | \
  awk '{{ if (substr($0,0,1) == "#") {{ print $0"\\toChrom\\toStart\\toEnd";}} else {{ print $0; }} }} ' | \
  {params.sd}/SelectBestcall.py --infile /dev/stdin --header int_filt_pb  --svlen 7  --key > {output.integrated_pbfinal}
cat {input.integrated} | \
  awk '{{ if (substr($0,0,1) == "#") {{ print "#sChrom\\tsStart\\tsEnd\\t"$0;}} else {{ print $1"\\t"$2"\\t"$3"\\t"$0;}} }}' | \
  bedtools slop -header -i stdin -b {params.window} -g {params.ref}.fai |\
  bedtools intersect -a stdin -b {input.pbfinalquery} -header -loj | \
  awk '{{ if (substr($0,0,1) == "#") {{ print $0"\\toChrom\\toStart\\toEnd";}} else {{ print $0; }} }} ' | \
  {params.sd}/SelectBestcall.py --infile /dev/stdin --header int_filt_pb  --svlen 7  --minOverlap 0.01 --key > {output.integrated_pbfinal}.any        
"""

#
#  Annotate pacbio calls with their overlap with Illumina calls.
#  Input (examples)
#     intfinal_query_pb_isect.DEL.bed - sloppy intersect, one line for each Illumina
#       call, with the key of the pacbio call if they overlap.
#     pbfinal.DEL.bed - the pacbio calls.
#  Output (examples)
#     pbfinal_int_tag.DEL.tab  - The key of the pacbio call this
#       Illumina call overlapped with.
#

rule AnnotatePBWithIntIsect:
    input:
        intfinalpbisect="intfinal_query_pb_isect.{svtype}.bed",
        pbfinal="pbfinal.{svtype}.bed"
    output:
        pbfinal_int_tag="pbfinal_int_tag.{svtype}.tab"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,
    shell:"""
bioawk -c hdr '{{ print $int_filt_pb_key;}}' < {input.intfinalpbisect}  | \
  {params.sd}/KeyToColumn.py --keys /dev/stdin --bed {input.pbfinal} --header int_match --out {output.pbfinal_int_tag}    
"""

rule CreateReverseLookupIntToPB:
    input:
        intfinalpbisect="intfinal_query_pb_isect.{svtype}.bed",
        pbfinal="pbfinal.{svtype}.bed"
    output:
        reverse="intfinal_reverse.{svtype}.tab"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,
    shell:"""
bioawk -c hdr '{{ print $int_filt_pb_key;}}' < {input.intfinalpbisect} | \
    {params.sd}/KeyToColumn.py --keys /dev/stdin --bed {input.pbfinal} --header int_match --out /dev/null --reverse {output.reverse}
"""



rule IntersectBNQueryAndIlluminaFiltered:
    input:
       intfilt="int_filt.{svtype}.bed",
       bnquery="bnquery.{svtype}.bed"
    output:
       intbn="int_filt_bn_isect.{svtype}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,
    shell:"""
header=`head -1 {input.bnquery}`
echo -e $header"\toChrom\toStart\toEnd" >  {output.intbn}.loj
bedtools intersect -a {input.bnquery} -b {input.intfilt} -loj >> {output.intbn}.loj
{params.sd}/SelectBestcall.py --infile {output.intbn}.loj --bnidx --ts 14 --te 15 --svlen 7 --out {output.intbn} --header int_isect --key
"""

#
# Given the full intersection of the pacbio final callset with the Illumina integrated,
# annotate which calls have the best overlap.
#
#rule BestPBFinalQueryAndIntFinal:
#    input:
#        pbfinalintisect="pbfinal_query_int_isect.{svtype}.bed"
#    output:
#        pbfinalintbest="pbfinal_query_int_best.{svtype}.bed"
#    params:
#        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
#        sd=SNAKEMAKE_DIR,
#    shell:"""
#header=`head -1 {input.pbfinalintisect}`
#echo -e $header"\toChrom\toStart\toEnd" >  {output.pbfinalintbest}
#{params.sd}/SelectBestcall.py --ts 35 --te 36 --infile {input.pbfinalintisect} --header intbest --key >> {output.pbfinalintbest}
#"""
#    

    
rule MakeSoftClip:
    input:
        bam=lambda wildcards:  localBams[wildcards.hap]
    output:
        clipped="clipped_hap.pb.h{hap}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,
        softclip=config["softclip"]
    shell:"""
if [ {params.softclip} == "NONE" ]; then
{params.sd}/PrintSoftClipped.py --bam {input.bam} --minClip 1000 --maxAlignLength 40000 --out {output.clipped}
else
cp {params.softclip}.{wildcards.hap}.bed {output.softclip}
"""

rule MakeILLSVCoverage:
    input:
        calls="integrated.{svtype}.bed"
    output:
        svCoverage="int_read_cov.{svtype}.bed",
    params:
        sge_opts="-pe serial 8 -l h_rt=1:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,
        fofn=config["fofn"],
        read_cov=config["bn_read_cov"]
    shell:"""
if [ {params.read_cov} == "NONE" ]; then
{params.sd}/../../utils/SVCoverage.py --fofn {params.fofn} --calls {input.calls} --out {output.svCoverage} --nproc 1 --bionano --op {wildcards.svtype} --header --window 1000
else
cp {params.read_cov}.{wildcards.svtype}.bed {output.svCoverage}
fi
"""

    
rule MakeBioNanoSVCoverage:
    input:
        bncalls="bncalls.{svtype}.bed"
    output:
        svCoverage="bn_read_cov.{svtype}.bed",
    params:
        sge_opts="-pe serial 8 -l h_rt=1:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,
        fofn=config["fofn"],
        read_cov=config["bn_read_cov"]
    shell:"""
if [ {params.read_cov} == "NONE" ]; then
{params.sd}/SVCoverage.py --fofn {params.fofn} --calls {input.bncalls} --out {output.svCoverage} --nproc 1 --bionano --op {wildcards.svtype} --header
else
cp {params.read_cov}.{wildcards.svtype}.bed {output.svCoverage}
fi
"""

        

rule CombineSoftClip:
    input:
        clippedAll=expand("clipped_hap.pb.h{hap}.bed",hap=haps),
    output:
        clippedCombined="clipped_all.pb.bed",
    params:
        sge_opts="-pe serial 8 -l h_rt=1:00:00 -l mfree=2G",
    shell:"""
cat {input.clippedAll} | bedtools sort > {output.clippedCombined}
"""
    
    
rule MakeReadIsect:
    input:
        sv="integrated.{svtype}.bed",
    output:
        readIsect="read_overlap.{svtype}.tab"
    params:
        gaps=config["pb_read_gaps"],
        sge_opts="-pe serial 8 -l h_rt=8:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,
        svOp=lambda wildcards: pbSVTypeMap[wildcards.svtype]
    shell:"""
cat {params.gaps} | bioawk -v op={params.svOp} -c hdr '{{ if (NR==1 || $svType == op) {{ print;}} }}' | \
bedtools intersect -sorted -f 0.5 -r -loj -a {input.sv} -b stdin | \
bedtools groupby  -g 1-5 -c 3 -o count | awk '{{ if (NR==1) print "number"; print $NF;}}' > {output.readIsect}
"""
rule MakeBNReadIsect:
    input:
        sv="bnquery.{svtype}.bed",
    output:
        readIsect="bn_read_overlap.{svtype}.bed"
    params:
        fofn=config["fofn"],
        sge_opts="-pe serial 8 -l h_rt=1:00:00 -l mfree=2G",
        sd=SNAKEMAKE_DIR,
        pbop=lambda wildcards: pbSVTypeMap[wildcards.svtype],
        bn_read_overlap=config["bn_read_overlap"]
    shell:"""
if [ {params.bn_read_overlap} == "NONE" ]; then
{params.sd}/FindRawReadSupport.py --calls {input.sv} --fofn {params.fofn} --out {output.readIsect} --nproc 1 --op {params.pbop} --isbn --header --ts 17 --te 18
else
cp {params.bn_read_overlap}.{wildcards.svtype}.bed {output.readIsect}
fi
"""

    

rule MakeBPQuery:
    input:
        sv="integrated.{svtype}.bed"
    output:
        query="bpquery.{svtype}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        ref=config["ref"],
        window=config["window"]
    shell:"""
cat {input.sv} | awk '{{ print $1"\\t"$2"\\t"$3"\\t"$0; }}' | bedtools slop -i stdin -b {params.window} -g {params.ref}.fai | awk '{{ print $1"\\t"$2"\\t"$2+2*{params.window}"\\t"$0; print $1"\\t"$3-2*{params.window}"\\t"$3"\\t"$0; }}' >  {output.query}
"""

rule MakeBNQuery:
    input:
        sv="bncalls.{svtype}.bed"
    output:
        query="bnquery.{svtype}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        ref=config["ref"],
        window=config["window"]
    shell:"""
nc=`head -1 {input.sv} | awk '{{ print NF+3;}}'`
cat {input.sv} | awk '{{ if (substr($0,0,1) == "#") {{ print "#sChrom\tsStart\tsEnd\t"$0; }} else {{ print $1"\t"$2"\t"$3"\t"$0;}} }}' | bedtools  slop -header -i stdin -b {params.window} -g {params.ref}.fai | bedtools groupby -header  -g 1-3 -c 4 -o first -full | cut -f 1-$nc > {output.query}
"""


        
        
rule MakeSVQuery:
    input:
        sv="integrated.{svtype}.bed"
    output:
        query="query.{svtype}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        ref=config["ref"],
        window=config["window"]
    shell:"""
cat {input.sv} | awk '{{ if (substr($0,0,1) == "#") {{ print "#qChrom\\tqStart\\tqEnd\\t"$0;}} else {{ print $1"\\t"$2"\\t"$3"\\t"$0; }} }}' | bedtools slop -header -i stdin -b {params.window} -g {params.ref}.fai > {output.query}
"""

rule OverlapBioNano:
    input:
        svQuery="query.{svtype}.bed",
        bnCalls="bncalls.{svtype}.bed"
    output:
        svIsect="bn_int_query_overlap.{svtype}.bed",
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        ref=config["ref"],
    shell:"""
header=`head -1 {input.svQuery}`
echo -e $header"\toChrom\toStart\toEnd" > {output.svIsect}
bedtools intersect -a {input.svQuery} -b {input.bnCalls} -loj >> {output.svIsect}
"""

rule SelectBestBNCall:
    input:
        svIsect="bn_int_query_overlap.{svtype}.bed",
    output:
        svBestCall="bnbest.{svtype}.tab",
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        sd=SNAKEMAKE_DIR
    shell:"""
{params.sd}/SelectBestcall.py --infile {input.svIsect} --header bnoverlap --bntarget > {output.svBestCall}
"""
        
rule OverlapIntegratedAndPB:
    input:
        svQuery="query.{svtype}.bed",
        pbPanel="pbcalls.{svtype}.bed",
    output:
        svIsect="query_overlap.{svtype}.bed",
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        ref=config["ref"],
    shell:"""
bedtools intersect -header -a  {input.svQuery} -b {input.pbPanel} -loj | awk '{{ if (substr($0,0,1) == "#") {{ print $0"\\toChrom\\toStart\\toEnd"; }} else {{ print $0; }} }}' > {output.svIsect}
"""
rule OverlapBNAndPB:
    input:
        svQuery="bnquery.{svtype}.bed",
        pbPanel="pbcalls.{svtype}.bed",
    output:
        svIsect="bnquery_ovp.{svtype}.bed",
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        ref=config["ref"],
    shell:"""
bedtools intersect -header -a {input.svQuery} -b {input.pbPanel} -loj > {output.svIsect}
"""

rule SelectBestPBCall:
    input:
        svIsect="query_overlap.{svtype}.bed",
    output:
        svBestCall="pbbest.{svtype}.tab",
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        sd=SNAKEMAKE_DIR
    shell:"""
{params.sd}/SelectBestcall.py --infile {input.svIsect} --header pbbest > {output.svBestCall}
"""
rule SelectBestPBCallForBN:
    input:
        svIsect="bnquery_ovp.{svtype}.bed",
    output:
        svBestCall="bnquery_best.{svtype}.tab",
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        sd=SNAKEMAKE_DIR
    shell:"""
cat {input.svIsect} | awk '{{ if (substr($0,0,1) == "#") {{ print $0"\toChrom\toStart\toEnd";}} else {{ print $0; }} }}' | \
 {params.sd}/SelectBestcall.py --infile /dev/stdin --bnidx --header pbbest >  {output.svBestCall} 
"""


rule CollectPBInclusive:
    input:
        localHap0=config["local"][0],
        localHap1=config["local"][1],
        stitchHap0=config["stitch"][0],
        stitchHap1=config["stitch"][1],
    output:
        mergedCalls="calls.pb.bed",
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
    shell:"""
head -1 {input.localHap0} | bioawk -c hdr '{{ print $_chrom"\\t"$tStart"\\t"$tEnd"\\t"$svType"\\t"$svLen"\\t"$qName"\\t"$qStart"\\t"$qEnd"\\tSOURCE\\tHAP";}}' > {output.mergedCalls}
cat {input.localHap0} | bioawk -c hdr '{{ print $_chrom"\\t"$tStart"\\t"$tEnd"\\t"$svType"\\t"$svLen"\\t"$qName"\\t"$qStart"\\t"$qEnd"\\tLOCAL\\tHAP0";}}' | grep -v "^#" >> {output.mergedCalls}
cat {input.localHap1} | bioawk -c hdr '{{ print $_chrom"\\t"$tStart"\\t"$tEnd"\\t"$svType"\\t"$svLen"\\t"$qName"\\t"$qStart"\\t"$qEnd"\\tLOCAL\\tHAP1";}}' | grep -v "^#" >> {output.mergedCalls}
cat {input.stitchHap0} | bioawk -c hdr '{{ print $_chrom"\\t"$tStart"\\t"$tEnd"\\t"$svType"\\t"$svLen"\\t"$qName"\\t"$qStart"\\t"$qEnd"\\tSTITCH\\tHAP0";}}' | grep -v "^#" >> {output.mergedCalls}
cat {input.stitchHap1} | bioawk -c hdr '{{ print $_chrom"\\t"$tStart"\\t"$tEnd"\\t"$svType"\\t"$svLen"\\t"$qName"\\t"$qStart"\\t"$qEnd"\\tSTITCH\\tHAP1";}}' | grep -v "^#" >> {output.mergedCalls}
"""

rule SeparatePBCallsByType:
    input:
        pbsv="calls.pb.bed",
    output:
        pbsvtype="pbcalls.{svtype}.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        pbop=lambda wildcards: pbSVTypeMap[wildcards.svtype]
    shell:"""
egrep "^#|{params.pbop}" {input.pbsv} | bedtools sort -header > {output.pbsvtype}
"""

    
rule MakeSVBed:
    input:
        svCallsVCF=config["svvcf"]
    output:
        svCallsBed="integrated.bed"
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        sd=SNAKEMAKE_DIR,
    shell:"""

{params.sd}/../../utils/variants_vcf_to_bed.py --vcf {input.svCallsVCF} --out {output.svCallsBed}.full --ignore-seqlen --filter --groupCommas  --fields CIPOS CIEND END SV_TYPE ALGORITHM --usesvlen
bioawk -c hdr '{{ print $_chrom"\\t"$tStart"\\t"$tEnd"\\t"$svType"\\t"$svLen"\\t"$svSeq"\\t"$CIEND"\\t"$CIPOS"\\t"$END"\\t"$SV_TYPE"\\t"$ALGORITHM"\\t"$FILTER }}' < {output.svCallsBed}.full > {output.svCallsBed}
"""

rule RenameIntegration:
    input:
        intg="integrated.bed"
    output:
        renamed="integrated_renamed.bed",
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        sd=SNAKEMAKE_DIR,        
    shell:"""
cat {input.intg} | {params.sd}/RenameEntries.py > {output.renamed}
"""

    
rule ExtractSVs:
    input:
        svcalls="integrated_renamed.bed"
    output:
        sepSV="integrated.{svtype}.bed",
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        sd=SNAKEMAKE_DIR,
    shell:"""
egrep "^#|{wildcards.svtype}" {input.svcalls} | bioawk  -c hdr '{{ if (NR==1 || ($3-$2 >= 50) ) print; }}' | cut -f 1-12 > {output.sepSV}
"""

rule NonRedundantToVCF:
    input:
        bed="{sample}.merged_nonredundant.{svtype}.bed",
    output:
        vcf="{sample}.merged_nonredundant.{svtype}.vcf",
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        sd=SNAKEMAKE_DIR,
        ref=config["ref"]
    shell:"""

{params.sd}/../../utils/variants_bed_to_vcf.py --bed {input.bed} --vcf {output.vcf} --type sv --ref {params.ref} --fields SVCLASS svClass CALLSET callset UNION union --sample {wildcards.sample}
"""

rule NonRedundantLociToVCF:
    input:
        bed="{sample}.merged_nonredundant_loci.DEL.bed",
    output:
        vcf="{sample}.merged_nonredundant_loci.DEL.vcf",
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        sd=SNAKEMAKE_DIR,
        ref=config["ref"]
    shell:"""
{params.sd}/../../utils/variants_bed_to_vcf.py --bed {input.bed} --vcf {output.vcf} --type sv --ref {params.ref} --fields SVCLASS svClass CALLSET callset UNION union --sample {wildcards.sample}
"""

rule RunVEP:
    input:
        vcf="{sample}.merged_nonredundant.{svtype}.vcf",
    output:
        vep="{sample}.merged_nonredundant_vep.{svtype}.vcf",
    params:
        sge_opts="-pe serial 6 -l h_rt=8:00:00 -l mfree=10G",
        sd=SNAKEMAKE_DIR,
        ref=config["ref"],
    shell:"""
{params.sd}/../../vep/RunVEP.sh {input.vcf} {output.vep}
"""
rule MakeVepTab:
    input:
        vep="{sample}.merged_nonredundant_vep.{svtype}.vcf",
    output:
        tab="{sample}.merged_nonredundant_vep.{svtype}.tab",
    params:
        sge_opts="-pe serial 6 -l h_rt=8:00:00 -l mfree=10G",
        sd=SNAKEMAKE_DIR,
        ref=config["ref"],
    shell:"""
{params.sd}/../../vep/VCF_VEP2TAB.pl {input.vep} > {output.tab}
"""
    
rule RunExons:
    input:
        bed="{sample}.merged_nonredundant.{svtype}.bed",
    output:
        vep="{sample}.merged_nonredundant_exon.{svtype}.bed",
    params:
        sge_opts="-pe serial 6 -l h_rt=2:00:00 -l mfree=8G",
        sd=SNAKEMAKE_DIR,
        sample=config["sample"],
        ref=config["ref"],
    shell:"""
{params.sd}/../../vep/intersect_exons.sh {input.bed} | bioawk -c hdr '{{ if (NR==1) {{ print $0"\\tsample";}}  else {{ print $0"\\t{params.sample}" }} }}' > {output.vep}
"""

rule AddPLI:
    input:
       bed="{sample}.merged_nonredundant_exon.{svtype}.bed",
    output:
       pli="{sample}.merged_nonredundant_exon.{svtype}.bed.pli",
    params:
        sge_opts="-pe serial 6 -l h_rt=2:00:00 -l mfree=8G",
        sd=SNAKEMAKE_DIR,
        sample=config["sample"],
        pli=config["pli"]
    shell:"""
{params.sd}/AddPLI.py --pli {params.pli} --svBed {input.bed} --out {output.pli}    
"""

    
rule AddRVIS:
    input:
        bed="{sample}.merged_nonredundant_exon.{svtype}.bed",
    output:
        rvis="{sample}.merged_nonredundant_exon.{svtype}.bed.rvis",
        rvistab="{sample}.merged_nonredundant_exon.{svtype}.bed.rvis.tab",        
    params:
        sge_opts="-pe serial 6 -l h_rt=2:00:00 -l mfree=8G",
        sd=SNAKEMAKE_DIR,
        sample=config["sample"],
        ref=config["ref"],
        rvis=config["rvis"]
    shell:"""
{params.sd}/AddRVIS.py --rvis {params.rvis} --svBed {input.bed} --out {output.rvis}.tab
paste {input.bed} {output.rvis}.tab > {output.rvis}
"""

rule MakeFullExonTable:
    input:
        pli="{sample}.merged_nonredundant_exon_pli.{svtype}.tab",    
        rvis="{sample}.merged_nonredundant_exon.{svtype}.bed.rvis.tab",
        bed="{sample}.merged_nonredundant_exon.{svtype}.bed",
        tr="{sample}.merged_nonredundant_exon.{svtype}.bed.tr",
    output:
        merged="{sample}.merged_nonredundant_exon.{svtype}.bed.merged",
    params:
        sge_opts="-pe serial 6 -l h_rt=2:00:00 -l mfree=8G",
        sd=SNAKEMAKE_DIR,
    shell:"""
paste {input.bed} {input.pli} {input.rvis} {input.tr} > {output.merged}
"""        

rule GetTFBSAblation:
    input:
        vep=config["sample"]+".merged_nonredundant_vep.DEL.vcf"
    output:
        tfbs="tfbs_ablation.bed"
    params:
        sd=SNAKEMAKE_DIR,
    shell:"""
{params.sd}/../..//vep/VCF_VEP2TAB.pl {input.vep} | egrep "^#|TFBS_ablati" | bedtools groupby -header -g  1-3 -c 4 -o first -full  > {output.tfbs}
"""
        
        
rule MakeCombinedVCF:
    input:
        mnr=[config["sample"]+".merged_nonredundant_loci.DEL.vcf", config["sample"]+".merged_nonredundant.INS.vcf"]
    output:
        comb=config["sample"] + ".merged_nonredundant.vcf"
    shell:"""
cat <( grep "^#" {input.mnr[0]} ) <( grep -v "^#" {input.mnr[0]} | sort -k1,1 -k2,2n ) <( grep -v "^#" {input.mnr[1]} | sort -k1,1 -k2,2n ) > {output.comb}
"""



