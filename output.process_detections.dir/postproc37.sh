#!/bin/bash
#collect catalogs obtained using different templates to provide a daily catalog
for day in 090330
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
