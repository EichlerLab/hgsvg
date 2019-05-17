#!/usr/bin/env bash

set -euo pipefail

bioawk -c hdr '{ if ($svType == "insertion" && NR > 1) { $tEnd=$tStart+$svLen;} print;}'
