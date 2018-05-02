from functions import detection as ds
infolder = 'Detections/'
alldet = ds.read_detections_cat(infolder)

for d in alldet:
    print(d.orig_time)
