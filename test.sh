chmod +x ./*.dir/*.py
cd ./download_data.dir
./download_eida_ingv.py
./download_inventory_ingv.py
cd ../trim_templates.dir
./trim_templates4.1.py
cd ../calc_ttimes.dir
./calcTT06.py
cd ../pympa38.dir
./pympa40mac.py 
cd ../postproc.dir
./postproc37.sh  
cd ../verify_detection.dir
./verify_detection03.py
