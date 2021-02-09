#!/bin/bash

# install GNU date "gdate" before running the script
# e.g. by calling "brew install coreutils" or "conda install coreutils"

# set start abd stop dates
start="2009-03-29"
stop="2009-03-30"

# generate a list of days to process
until [[ $start > $stop ]]; do
    echo $start
    list+=`echo "${start:2:2}${start:5:2}${start:8:2} "`
    start=$(date -I -d "$start+ 1 day")
done

#collect catalogs obtained using different templates to provide a daily catalog
for day in $list
do
find ./ -maxdepth 1 -name "*.${day}.cat" -print0 | xargs -0 cat > ${day}Acat

# filter and cut daily catalog to select maximum threshold events
#./filterCATtimestamp.py<<END
rm dcat
cp ${day}Acat dcat
./process_detections.py
cp dcatf1f2 ${day}Acatf1f2
rm dcatf1 dcatf1f2
done

# get <<outcat>> the  cleaned catalogue
cat *Acatf1f2 > outcat
