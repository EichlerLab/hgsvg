#!/usr/bin/env python



import argparse
import os
import subprocess
import sys

ap = argparse.ArgumentParser(description="Print the regions to map")
ap.add_argument("--rgn", help="Regions file", required=True)
ap.add_argument("--params",help="Parameter file.",required=True)
ap.add_argument("--regionsPerJob",help="Run this many regions in each job.", default=15,type=int)
ap.add_argument("--assemblies", help="Assemblies directory",default="assemblies")
ap.add_argument("--wd", help="Working directory, default=cwd",default=None)
ap.add_argument("--tc", help="TC (number of concurrent jobs", default=300)
ap.add_argument("--mem", help="Memory per job", default="4")
ap.add_argument("--trio",help="Run on trio instead of single sample", action='store_true', default=False)
ap.add_argument("--joblimit",help="Maximum number of jobs to have in one run.",type=int,default=75000)
ap.add_argument("--noDelay", help="Immediately start jobs rather than staggered start", action='store_true', default=False)
ap.add_argument("--slots", help="Number of slots to reserve", default=4,type=int)


if os.path.exists("config.json") == False:
    print("ERROR. The file 'config.json' for tiled-meta assembly must exist.")
    sys.exit(1)

args = ap.parse_args()

rgnFile = open(args.rgn)
regions = [l.rstrip() for l in rgnFile]

jobNumber = 0
rgnNumber = 0
regionsFileName = "regions.txt"
regionsFile = open(regionsFileName,'w')
if args.wd is None:
    args.wd = os.getcwd()

args.driverScript="RunTiledAssemblyOnRegions.sh"
args.makefile="RunPhasePartitionedAssemblies.mak"
if args.trio is True:
    args.driverScript="RunTrioTiledAssemblyOnRegions.sh"
    args.makefile="RunTrioPartionedAssemblies.mak"

cmdStr ="""#!/usr/bin/env  bash
#
#$ -t 1-{} -tc {} -S /bin/bash -V  -e /dev/null -o /dev/null -l mfree={}G -l h_rt=02:00:00 -pe serial {} -wd {} -p -600 -q eichler-short.q
#
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/local_assembly/{} `awk "NR == $SGE_TASK_ID" {} ` lasm \
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/local_assembly/{}  {} $SGE_TASK_ID {}
"""    

def PrintJobsCommand(outFile, nJobs, rgnFile, args):
    outFile.write(cmdStr.format(nJobs, args.tc, args.mem, args.slots, args.wd, args.driverScript, rgnFile, args.makefile, args.params,args.delay))
    outFile.close()


remainingRegions = []
assemblyFofn = open("assemblies.fofn",'w')
cmd = "mkdir -p "  + args.assemblies
subprocess.call(cmd.split())
cmd = "ls {}/".format(args.assemblies)
p=subprocess.Popen(cmd.split(),stdout=assemblyFofn)
p.wait()
assemblyFofn.close()
assemblyFofn = open("assemblies.fofn")
assemblies = { l[0:l.rfind(".")] for l in assemblyFofn }

for r in regions:
    if r not in assemblies:
        if "chrY" not in r:
            remainingRegions.append(r)

i = 0

def GetChrom(rgn):
    dot = rgn.find('.')
    return rgn[0:dot]

nRegions = 0
while i < len(remainingRegions):
    j=i
    while j < len(remainingRegions) and j < i+args.regionsPerJob and GetChrom(remainingRegions[i]) == GetChrom(remainingRegions[j]):
        j+=1
    rgnGroup = ";".join(remainingRegions[i:j])
    regionsFile.write(rgnGroup + "\n")
    nRegions+=1
    i=j

cmdFile = open("AssembleRegions.sh",'w')
if args.noDelay:
    args.delay = 0
else:
    args.delay = args.tc


PrintJobsCommand(cmdFile, nRegions, regionsFileName, args)
cmdFile.close()

submit = open("submit_jobs.sh",'w')
jobIndex = 0
submit.write("res0=`~mchaisso/scripts/qsubid -cwd ./AssembleRegions.sh`\n")
submit.write("res1=`~mchaisso/scripts/qsubid -hold_jid $res0 -cwd ./AssembleRegions.sh`\n")
submit.write("res2=`~mchaisso/scripts/qsubid -hold_jid $res1 -cwd ./AssembleRegions.sh`\n")
submit.write("res1=`~mchaisso/scripts/qsubid -hold_jid $res2 -cwd -l h_rt=24:00:00 -pe serial 8 -l mfree=1G /net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg/sv/run_diploid_annotation.sh {}`\n".format(args.params))
submit.write("res2=`~mchaisso/scripts/qsubid -hold_jid $res3 -cwd -l h_rt=48:00:00 -pe serial 1 -l mfree=1G /net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg/stitching/run_snakemake.sh`\n")
submit.close()
