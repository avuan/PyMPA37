class detection_cat:
    def __init__(self, temp_used, orig_time, mag, avr_cc, thresh_val, avr_cc_no_shift, thresh_val_no_shift, no_chan):
        self.temp_used = temp_used
        self.orig_time = orig_time
        self.mag = mag
        self.avr_cc = avr_cc
        self.thresh_val = thresh_val
        self.avr_cc_no_shift = avr_cc_no_shift
        self.thresh_val_no_shift = thresh_val_no_shift
        self.no_chan = no_chan

def read_detections_cat(infolder):
    import csv
    from obspy import UTCDateTime
    import os
    
    filesfolder = os.listdir(infolder)
    all_det = []
    
    for detfile in filesfolder:
        if detfile.endswith('.cat'):
            with open(infolder+detfile) as detections_file:
                detections = csv.reader(detections_file, delimiter=' ')
                
                for line in detections:
                    temp_used = int(line[0])
                    orig_time = UTCDateTime(line[1])
                    mag = float(line[2])
                    avr_cc = float(line[3])
                    thresh_val = float(line[4])
                    avr_cc_no_shift = float(line[5])
                    thresh_val_no_shift = float(line[6])
                    no_chan = float(line[7])
                    
                    det = detection_cat(temp_used,
                                        orig_time,
                                        mag,
                                        avr_cc,
                                        thresh_val,
                                        avr_cc_no_shift,
                                        thresh_val_no_shift,
                                        no_chan)
                    all_det.append(det)
    
                detections_file.close()
    return all_det

class detection_stats:
    def __init__(self, station, channel, xcorr, xcorr_no_shift, shift_in_samples, day, temp_used, orig_time, mag_det, mag_temp, no_chan_above_xcorr_35, mad, avr_cc, thresh_val, avr_cc_no_shift, thresh_val_no_shift, no_chan_above_xcorr_30, no_chan_above_xcorr_50, no_chan_above_xcorr_70, no_chan_above_xcorr_90):
        
        self.station = station
        self.channel = channel
        self.xcorr = xcorr
        self.xcorr_no_shift = xcorr_no_shift
        self.shift_in_samples = shift_in_samples
        self.day = day
        self.temp_used = temp_used
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
    all_det = []
    
    for detfile in filesfolder:
        if detfile.endswith('.stats'):
            with open(infolder+detfile) as detections_file:
                detections = detections_file.read()
                detections = detections.split(';')
                for detection in detections:
                    line = detection.split('\n')
                    line = list(filter(None, line))
                    
                    station = None
                    channel = None
                    xcorr = None
                    xcorr_no_shift = None
                    shift_in_samples = None
                    day = None
                    temp_used = None
                    orig_time = None
                    mag_det = None
                    mag_temp = None
                    no_chan_above_xcorr_35 = None
                    mad = None
                    avr_cc = None
                    thresh_val = None
                    avr_cc_no_shift = None
                    thresh_val_no_shift = None
                    no_chan_above_xcorr_30 = None
                    no_chan_above_xcorr_50 = None
                    no_chan_above_xcorr_70 = None
                    no_chan_above_xcorr_90 = None
                    
                    for element in line:
                        if not element[0:2].isdigit():
                            ele = element.split(' ')
                            station = str(ele[0])
                            channel = str(ele[1])
                            xcorr = float(ele[2])
                            xcorr_no_shift = float(ele[3])
                            shift_in_samples = float(ele[4])
                            print(station, channel, xcorr, xcorr_no_shift, shift_in_samples)
                        else:
                            ele = element.split(' ')
                            day = int(ele[0])
                            temp_used = int(ele[1])
                            unknown = int(ele[2])
                            orig_time = UTCDateTime(ele[3])
                            mag_det = float(ele[4])
                            mag_temp = float(ele[5])
                            no_chan_above_xcorr_35 = float(ele[6])
                            mad = float(ele[7])
                            avr_cc = float(ele[8])
                            thresh_val = float(ele[9])
                            avr_cc_no_shift = float(ele[10])
                            thresh_val_no_shift = float(ele[11])
                            no_chan_above_xcorr_30 = float(ele[12])
                            no_chan_above_xcorr_50 = float(ele[13])
                            no_chan_above_xcorr_70 = float(ele[14])
                            no_chan_above_xcorr_90 = float(ele[15])
                
                    print(day, temp_used, unknown, orig_time, mag_det, mag_temp, no_chan_above_xcorr_35, mad, avr_cc, thresh_val, avr_cc_no_shift, thresh_val_no_shift, no_chan_above_xcorr_30, no_chan_above_xcorr_50, no_chan_above_xcorr_70, no_chan_above_xcorr_90, '\n')
                    
                    det = detection_stats(station,
                                          channel,
                                          xcorr,
                                          xcorr_no_shift,
                                          shift_in_samples,
                                          day,
                                          temp_used,
                                          orig_time,
                                          mag_det,
                                          mag_temp,
                                          no_chan_above_xcorr_35,
                                          mad,
                                          avr_cc,
                                          thresh_val,
                                          avr_cc_no_shift,
                                          thresh_val_no_shift,
                                          no_chan_above_xcorr_30,
                                          no_chan_above_xcorr_50,
                                          no_chan_above_xcorr_70,
                                          no_chan_above_xcorr_90)
                    all_det.append(det)
        
            detections_file.close()
    return all_det

