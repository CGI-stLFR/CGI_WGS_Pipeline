#!/usr/bin/env bash

awk -F $'\\t' '(NR>1){fragcount++;fraglength=fraglength+$6;fragreadcount=fragreadcount+$7}END{ave_fraglen=fraglength/fragcount;ave_fragreadcount=fragreadcount/fragcount;print "fragcount","ave_fraglen","ave_fragreadcount";print fragcount,ave_fraglen,ave_fragreadcount}' OFS="\t" $1 > $2
