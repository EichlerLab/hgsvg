#!/usr/bin/env bash

set -euo pipefail

c=0
IFS=$'\n'

while read line; do
		tmp=`echo $line | sed 's/\\t/\\\\t/'`"\n";
		lines[$c]=$tmp
		c=$(expr $c + 1);
done
if [ ${#lines[@]} -gt 1 ]; then 
		echo -e ${lines[*]} | sed 's/^ //g' |  sed '/^$/d'  | eval $1
else
		echo -e ${lines[*]} | sed 's/^ //g' |  sed '/^$/d' 
fi
