class detection:
    def __init__(self, temp_used, orig_time, mag, avr_cc, thresh_val, avr_cc_no_shift, thresh_val_no_shift, no_chan):
        self.temp_used = temp_used
        self.orig_time = orig_time
        self.mag = mag
        self.avr_cc = avr_cc
        self.thresh_val = thresh_val
        self.avr_cc_no_shift = avr_cc_no_shift
        self.thresh_val_no_shift = thresh_val_no_shift
        self.no_chan = no_chan

def read_detections(infile):
    import csv
    from obspy import UTCDateTime
    all_det = []
    
    with open(infile) as detections_file:
        detections = csv.reader(detections_file, delimiter='\t')
        headers = next(detections)
        
        for line in detections:
            temp_used = int(line[0])
            orig_time = UTCDateTime(line[1])
            mag = float(line[2])
            avr_cc = float(line[3])
            thresh_val = float(line[4])
            avr_cc_no_shift = float(line[5])
            thresh_val_no_shift = float(line[6])
            no_chan = float(line[7])
            
            det = detection(temp_used,
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

def decluster(all_det, between_time):
    unique_detections = []
    for master in all_det:
        keep = True
        for slave in all_det:
            if not master == slave and abs(master.orig_time - slave.orig_time) <= float(between_time):
                if not master.thresh_val >= slave.thresh_val:
                    keep = False
                    break
        if keep:
            unique_detections.append(master)
    return unique_detections

def mag_time_plot(detections):
    from matplotlib import pyplot as plt
    mags = []
    times = []
    for detection in detections:
        mags.append(detection.mag)
        times.append(detection.orig_time.datetime)
    minmag = min(mags)
    maxmag = max(mags)

    plt.plot(times, mags, 'ko')
    plt.ylim(minmag-0.5, maxmag+0.5)
    plt.show()

def high_cc(all_det, ccval):
    high_cc_val = []
    for detection in all_det:
        if detection.avr_cc >= float(ccval):
            high_cc_val.append(detection)
    else:
        pass
    return high_cc_val

def extract_detections(detections, st_path, freqmin, freqmax):
    from obspy import read
    for detection in detections:
        det_year = str(detection.orig_time.year)[2:]
        det_month = str(detection.orig_time.month)
        det_day = str(detection.orig_time.day)
        
        if int(det_month) < 10:
            det_month = '0' + det_month
        else:
            det_month = det_month
        
        if int(det_day) < 10:
            det_day = '0' + det_day
        else:
            det_day = det_day

        str_time = det_year + det_month + det_day
        
        st = read(st_path+str_time+'.*')
        st.merge()
        print(st)
        st.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
        st.trim(starttime=detection.orig_time+2.0, endtime=detection.orig_time+15.0)
        st.write(str(detection.orig_time)+'.mseed', format='MSEED')
        st.normalize()
        st.plot(outfile=str(detection.orig_time)+'.png', format='PNG')
        print('Done with :', str(detection.orig_time))


