#!/usr/bin/env python

import argparse

ap = argparse.ArgumentParser(description="Collect")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()
