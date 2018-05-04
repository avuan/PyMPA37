from functions import detection_stats as ds

infolder = 'Detections/'
alldets = ds.read_detections_stats(infolder)

for info in alldets:
    sta_chan_stats, origins = info
    for sta_chan_stat in sta_chan_stats:
        if sta_chan_stat.xcorr > 0.5:
            print(sta_chan_stat.station)

    print(origins.orig_time)
