#!/usr/bin/env bash

awk -F $'\t' '(NR>1){barcodecount++;fragonbarcode=fragonbarcode+$2;barcodereadcount=barcodereadcount+$3}END{ave_fragonbarcode=fragonbarcode/barcodecount;ave_barcodereadcount=barcodereadcount/barcodecount;print "barcodecount","ave_fragonbarcode","ave_barcodereadcount";print barcodecount,ave_fragonbarcode,ave_barcodereadcount}' OFS="\t" $1 > $2
