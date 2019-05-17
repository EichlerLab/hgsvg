#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -ne 2 ]; then
   echo "Usage: MakeManifest.sh file_or_directory manifest.tsv"
   exit 1
fi
echo $1 >> /dev/stderr
file-manifest $1 > $2
