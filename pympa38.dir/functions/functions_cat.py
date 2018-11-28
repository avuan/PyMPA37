def decluster_cat(all_det, between_time):
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

def mag_time_plot_cat(detections):
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

def high_cc_cat(all_det, ccval):
    high_cc_val = []
    for detection in all_det:
        if detection.avr_cc >= float(ccval):
            high_cc_val.append(detection)
    else:
        pass
    return high_cc_val

def extract_detections_cat(detections, st_path, freqmin, freqmax):
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
