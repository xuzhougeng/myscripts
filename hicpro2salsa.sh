#!/usr/bin/bash

input=$1
output=$2

awk 'BEGIN{OFS="\t"} {print $2,$3,$3+50,$1"/1","60",$4} {print $5,$6,$6+50,$1"/2","60",$7}' $1 > $2
