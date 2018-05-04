class sta_chan_stat:
    def __init__(self, station, channel, xcorr, xcorr_no_shift, shift_in_samples):

        self.station = station
        self.channel = channel
        self.xcorr = xcorr
        self.xcorr_no_shift = xcorr_no_shift
        self.shift_in_samples = shift_in_samples

class detection_stats:
    def __init__(self, day, temp_used, unknown, orig_time, mag_det, mag_temp, no_chan_above_xcorr_35, mad, avr_cc, thresh_val, avr_cc_no_shift, thresh_val_no_shift, no_chan_above_xcorr_30, no_chan_above_xcorr_50, no_chan_above_xcorr_70, no_chan_above_xcorr_90):

        self.day = day
        self.temp_used = temp_used
        self.unknown = unknown
        self.orig_time = orig_time
        self.mag_det = mag_det
        self.mag_temp = mag_temp
        self.no_chan_above_xcorr_35 = no_chan_above_xcorr_35
        self.mad = mad
        self.avr_cc = avr_cc
        self.thresh_val = thresh_val
        self.avr_cc_no_shift = avr_cc_no_shift
        self.thresh_val_no_shift = thresh_val_no_shift
        self.no_chan_above_xcorr_30 = no_chan_above_xcorr_30
        self.no_chan_above_xcorr_50 = no_chan_above_xcorr_50
        self.no_chan_above_xcorr_70 = no_chan_above_xcorr_70
        self.no_chan_above_xcorr_90 = no_chan_above_xcorr_90

def read_detections_stats(infolder):
    import csv
    from obspy import UTCDateTime
    import os
    
    filesfolder = os.listdir(infolder)
    
    for detfile in filesfolder:
        if detfile.endswith('.stats'):
            with open(infolder+detfile) as detections_file:
                detections = detections_file.read()
                detections = detections.split(';')
                
                all_dets = []
                
                for detection in detections:
                    
                    stachanstats = []

                    detection = detection.split('\n')
                    detection = list(filter(None, detection))

                    for line in detection:
                        if not (line[0:2]).isdigit():
                            
                            station_element = line.split(' ')
                            station_element = list(filter(None, station_element))
                            station = str(station_element[0])
                            channel = str(station_element[1])
                            xcorr = float(station_element[2])
                            xcorr_no_shift = float(station_element[3])
                            shift_in_samples = float(station_element[4])
                            
                            stachanstat = sta_chan_stat(station,
                                                        channel,
                                                        xcorr,
                                                        xcorr_no_shift,
                                                        shift_in_samples)
                            stachanstats.append(stachanstat)
                            
                        else:
                            origin_element = line.split(' ')
                            origin_element = list(filter(None, origin_element))
                            day = int(origin_element[0])
                            temp_used = int(origin_element[1])
                            unknown = int(origin_element[2])
                            orig_time = UTCDateTime(origin_element[3])
                            mag_det = float(origin_element[4])
                            mag_temp = float(origin_element[5])
                            no_chan_above_xcorr_35 = float(origin_element[6])
                            mad = float(origin_element[7])
                            avr_cc = float(origin_element[8])
                            thresh_val = float(origin_element[9])
                            avr_cc_no_shift = float(origin_element[10])
                            thresh_val_no_shift = float(origin_element[11])
                            no_chan_above_xcorr_30 = float(origin_element[12])
                            no_chan_above_xcorr_50 = float(origin_element[13])
                            no_chan_above_xcorr_70 = float(origin_element[14])
                            no_chan_above_xcorr_90 = float(origin_element[15])
                            
                            origin_data = detection_stats(day, temp_used, unknown, orig_time, mag_det, mag_temp, no_chan_above_xcorr_35, mad, avr_cc, thresh_val, avr_cc_no_shift, thresh_val_no_shift, no_chan_above_xcorr_30, no_chan_above_xcorr_50, no_chan_above_xcorr_70, no_chan_above_xcorr_90)

                    det_tuple = (stachanstats, origin_data)
                    all_dets.append(det_tuple)

                return(all_dets)

