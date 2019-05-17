#!/usr/bin/env bash

set -euo pipefail

bioawk -c hdr 'BEGIN{OFS="\t";} { if ($svType == "insertion" && NR > 1) { $tEnd=$tStart+1;} print;}'
