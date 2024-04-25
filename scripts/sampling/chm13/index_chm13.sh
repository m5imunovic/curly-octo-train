#!/bin/bash

mkdir paf
REF_DIR=../references/homo_sapiens/chm13_v2
minimap2 -x map-hifi -H -t 20 -d ./ref.mmi "$REF_DIR"/chm13v2.0.fa.gz
