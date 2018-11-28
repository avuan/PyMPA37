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
