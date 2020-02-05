#!/usr/bin/env bash

echo barcode$'\t'total_min5000frag_count$'\t'total_min5000frag_readcount > $2

awk '{if(NR>1){
barcode=$1
fragreadcount=$7

if(prev_barcode==""){
fraginbarcodecount=1
totalreadcount=fragreadcount}

else if(barcode==prev_barcode){
fraginbarcodecount++
totalreadcount=totalreadcount+fragreadcount}

else if(barcode!=prev_barcode){
print prev_barcode,fraginbarcodecount,totalreadcount
fraginbarcodecount=1
totalreadcount=fragreadcount
}

prev_barcode=barcode

}}' OFS="\t" $1 >> $2 
