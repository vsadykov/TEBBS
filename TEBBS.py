import numpy
import datetime
import os, sys
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import requests
import wget
import shutil

# function to read the Temperature-Denomenator array from the file for the corresponding GOES satellite
def read_t_denom_file(gver):
    str_gver = str(gver)
    if (len(str_gver) < 2): str_gver = '0'+str_gver
    # opening temp-denom file
    fcur = open('temperature_denom_g'+str_gver+'.txt')
    # finding number of strings in the file
    t_denom = []
    for line in fcur:
        t = float(line[0:16])
        denom = float(line[16:-1])
        t_denom.append([t,denom])
    t_denom = numpy.array(t_denom)
    return t_denom


# function to read the Ratio-Temperature array from the file for the corresponding GOES satellite
def read_ratio_t_file(gver):
    str_gver = str(gver)
    if (len(str_gver) < 2): str_gver = '0'+str_gver
    # opening temp-denom file
    fcur = open('temperature_ratio_g'+str_gver+'.txt')
    # finding number of strings in the file
    ratio_t = []
    for line in fcur:
        ratio = float(line[0:16])
        t = float(line[16:-1])
        ratio_t.append([ratio,t])
    ratio_t = numpy.array(ratio_t)
    return ratio_t


# function to receive flux Ratio from the flux values
def get_ratio(Aflux,Bflux):
    ratio = 0.82353*Aflux/Bflux
    return ratio


# function to receive EM from the Bflux and Denomenator
def get_em(Bflux,denom):
    em = (1000000.0/0.7)*Bflux/denom
    return em


# function to convert datetime string to python datetime function
def return_datetime(inputtime):
    ctime = datetime.datetime.now()
    ctime = ctime.replace(year = int(inputtime[0:4]))
    ctime = ctime.replace(month = int(inputtime[5:7]))
    ctime = ctime.replace(day = int(inputtime[8:10]))
    ctime = ctime.replace(hour = int(inputtime[11:13]))
    ctime = ctime.replace(minute = int(inputtime[14:16]))
    ctime = ctime.replace(second = int(inputtime[17:19]))
    return ctime


# function to convert datetime string to seconds from the beginning of the day
def return_sec(inputtime):
    hour = int(inputtime[11:13])
    minute = int(inputtime[14:16])
    sec = int(inputtime[17:19])
    sec += minute*60 + hour*3600
    return sec


# function to read GOES fluxes using the input file
def read_flux(gfile):
    hdulist = fits.open(gfile)
    if (type(hdulist[0].data) is not numpy.ndarray):
        data = hdulist[2].data[0]
        timing = data[0]
        fluxes = data[1]
        return timing, fluxes
    else:
        timing = numpy.zeros([hdulist[0].data.shape[1]], dtype='>f8')
        timing = hdulist[0].data[0].astype('float64')
        fluxes = numpy.zeros([hdulist[0].data.shape[1],2], dtype='>f4')
        fluxes[:,0] = hdulist[0].data[1].astype('float32')
        fluxes[:,1] = hdulist[0].data[2].astype('float32')
        return timing, fluxes


# function to find the maximum time for the flare flux
def find_max_sec(timing, fluxes, flare_start_time, flare_end_time):
    maxtime_Aflux = 0
    maxtime_Bflux = 0
    max_Aflux = 0
    max_Bflux = 0
    for i in range (0,timing.shape[0],1):
        if ((timing[i] >= flare_start_time) and (timing[i] <= flare_end_time)):
            if (fluxes[i,0] > max_Bflux):
                max_Bflux = fluxes[i,0]
                maxtime_Bflux = timing[i]
            if (fluxes[i,1] > max_Aflux):
                max_Aflux = fluxes[i,1]
                maxtime_Aflux = timing[i]
    if (maxtime_Aflux > maxtime_Bflux):
        return maxtime_Bflux
    else:
        return maxtime_Aflux


# function to extract minimum values (lower bounds) for fluxes from start time till peak time
def extract_min(timing, flare_start_time, flare_peak_time, fluxes):
    start_index = numpy.argmin(abs(timing - flare_start_time))
    end_index = numpy.argmin(abs(timing - flare_peak_time))
    Amin = numpy.amin(fluxes[start_index:end_index,1])
    Bmin = numpy.amin(fluxes[start_index:end_index,0])
    return Amin, Bmin


# function to define background grid
def make_grid(Amin, Bmin):
    bgrid = numpy.zeros([20,20,2])
    for i in range (0,20,1):
        for j in range (0,20,1):
            bgrid[i,j,0] = float(i+0.5)*Amin/20.0
            bgrid[i,j,1] = float(j+0.5)*Bmin/20.0
    return bgrid


# function to extract backround-subtracted fluxes for the rising phase of the flare (without first 1/6 as suggested in the paper)
def extract_fluxes(timing, flare_start_time, flare_peak_time, fluxes, bgrid):
    start = flare_start_time + (flare_peak_time-flare_start_time)/6
    end = flare_peak_time
    start_index = numpy.argmin(abs(timing - start))
    start_index_ext = numpy.argmin(abs(timing - flare_start_time))
    end_index = numpy.argmin(abs(timing - end))
    # arrays with fluxes
    fluxes_bsub = numpy.zeros([end_index-start_index,bgrid.shape[0],bgrid.shape[1],bgrid.shape[2]], dtype = float)
    fluxes_bsub_ext = numpy.zeros([end_index-start_index_ext,bgrid.shape[0],bgrid.shape[1],bgrid.shape[2]], dtype = float)
    for i in range (0,bgrid.shape[0],1):
        for j in range (0,bgrid.shape[1],1):
            fluxes_bsub[:,i,j,0] = fluxes[start_index:end_index,1] - bgrid[i,j,0]
            fluxes_bsub[:,i,j,1] = fluxes[start_index:end_index,0] - bgrid[i,j,1]
            fluxes_bsub_ext[:,i,j,0] = fluxes[start_index_ext:end_index,1] - bgrid[i,j,0]
            fluxes_bsub_ext[:,i,j,1] = fluxes[start_index_ext:end_index,0] - bgrid[i,j,1]
    return fluxes_bsub, fluxes_bsub_ext, timing[start_index:end_index], timing[start_index_ext:end_index]


# function to calculate Temperature and EM for the fluxes array; the last dimension assume to have 2 values corr. to the fluxes
def calculate_temperature_em(fluxes, gnum):
    atemp = numpy.zeros([fluxes.shape[0],fluxes.shape[1],fluxes.shape[2]])
    aem = numpy.zeros([fluxes.shape[0],fluxes.shape[1],fluxes.shape[2]])
    ratio_t = read_ratio_t_file(gnum)
    t_denom = read_t_denom_file(gnum)
    #ftemp = interp1d(ratio_t[:,0],ratio_t[:,1],kind='quadratic',bounds_error=False)
    #fdenom = interp1d(t_denom[:,0],t_denom[:,1],kind='quadratic',bounds_error=False)
    for t in range (0,fluxes.shape[0]):
        for i in range (0,fluxes.shape[1]):
            for j in range (0,fluxes.shape[2]):
                Aflux = fluxes[t,i,j,0]
                Bflux = fluxes[t,i,j,1]
                ratio = get_ratio(Aflux,Bflux)
                temp = numpy.interp(ratio,ratio_t[:,0],ratio_t[:,1])
                denom = numpy.interp(temp,t_denom[:,0],t_denom[:,1])
                #if (ratio < ratio_t[0,0]): ratio = ratio_t[0,0]
                #if (ratio > ratio_t[-1,0]): ratio = ratio_t[-1,0]
                #temp = ftemp(ratio)
                #denom = fdenom(temp)
                em = get_em(Bflux,denom)
                atemp[t,i,j] = temp
                aem[t,i,j] = em
    return atemp, aem


# hot flare test; if the temperature is < 4MK for one of the points of the flare, discard the background combination
def test_temperature_min(temp, labelgrid):
    temperature_min_flag = 1
    labelgrid_copy = numpy.copy(labelgrid)
    for i in range (0,temp.shape[1],1):
        for j in range (0,temp.shape[2],1):
            if (numpy.amin(temp[:,i,j]) < 4.0): labelgrid[i,j] = 0
    # if there are no solutions satisfying this condition, let's leave all of them
    if (sum(sum(labelgrid)) == 0):
        labelgrid = numpy.copy(labelgrid_copy)
        temperature_min_flag = 0
    return labelgrid, temperature_min_flag


# increasing array test; if the array growth is < average -7%, discard the background combination
def test_array_increase(temp, labelgrid):
    total_increases = 0
    for i in range (0,temp.shape[1],1):
        for j in range (0,temp.shape[2],1):
            if (labelgrid[i,j] == 1):
                for t in range (1,temp.shape[0],1):
                    if (temp[t,i,j] > temp[t-1,i,j]): total_increases += 1
    total_increases = total_increases/sum(sum(labelgrid))   # normalizing to number of curves
    total_increases -= int(0.07*float(total_increases))    # -7% as suggested in the paper
    for i in range (0,temp.shape[1],1):
        for j in range (0,temp.shape[2],1):
            if (labelgrid[i,j] == 1):
                increases = 0
                for t in range (1,temp.shape[0],1):
                    if (temp[t,i,j] > temp[t-1,i,j]): increases += 1
                if (increases <= total_increases): labelgrid[i,j] = 0
    return labelgrid


# test for the initial 1/6 of the growing phase. If the peak in extended flux is greater than the peak in removed 1/5th one,
# discard the background combination. If the output is zero, return one with smallest peak deviation
def test_t_em_initials(temp, temp_ext, em, em_ext, labelgrid):
    initials_flag = 1
    labelgrid_temp = numpy.copy(labelgrid)
    labelgrid_em = numpy.copy(labelgrid)
    labelgrid_ext = numpy.copy(labelgrid)
    for i in range (0,temp.shape[1],1):
        for j in range (0,temp.shape[2],1):
            if (numpy.amax(temp_ext[:,i,j]) > numpy.amax(temp[:,i,j])): labelgrid_temp[i,j] = 0
            if (numpy.amax(em_ext[:,i,j]) > numpy.amax(em[:,i,j])): labelgrid_em[i,j] = 0
    # if there are no solutions in the parameter space, find one with the lowest deviation. Separated for T and EM.
    indexes = [-1,-1]
    for i in range (0,temp.shape[1],1):
            for j in range (0,temp.shape[2],1):
                labelgrid_ext[i,j] = labelgrid[i,j]*labelgrid_temp[i,j]*labelgrid_em[i,j]
    if (sum(sum(labelgrid_ext)) == 0):
        initials_flag = 0
        if ((sum(sum(labelgrid_temp)) == 0) and sum(sum(labelgrid_em)) > 0):
            for i in range (0,temp.shape[1],1):
                for j in range (0,temp.shape[2],1):
                    if (labelgrid[i,j] == 1): indexes = [i,j]
            deviation = numpy.amax(temp_ext[:,indexes[0],indexes[1]]) - numpy.amax(temp[:,indexes[0],indexes[1]])
            for i in range (0,temp.shape[1],1):
                for j in range (0,temp.shape[2],1):
                    dev = numpy.amax(temp_ext[:,i,j]) - numpy.amax(temp[:,i,j])
                    if ((dev < deviation) and (labelgrid[i,j] == 1)):
                        deviation = dev
                        indexes = [i,j]
        else:
            for i in range (0,em.shape[1],1):
                for j in range (0,em.shape[2],1):
                    if (labelgrid[i,j] == 1): indexes = [i,j]
            deviation = numpy.amax(em_ext[:,indexes[0],indexes[1]]) - numpy.amax(em[:,indexes[0],indexes[1]])
            for i in range (0,em.shape[1],1):
                for j in range (0,em.shape[2],1):
                    dev = numpy.amax(em_ext[:,i,j]) - numpy.amax(em[:,i,j])
                    if ((dev < deviation) and (labelgrid[i,j] == 1)):
                        deviation = dev
                        indexes = [i,j]
    # if everything is all right
    if (indexes != [-1,-1]):
        labelgrid_ext[indexes[0],indexes[1]] = 1
    return labelgrid_ext, initials_flag


# define error ranges for the selected array
def error_range(inputarray,labelgrid):
    cparray = numpy.copy(inputarray)
    # pushing first 1/6 and last 1/10 cparray values to zero
    for i in range (0,cparray.shape[0]/6,1):
        cparray[i,:,:] *= 0.0
    for i in range (9*cparray.shape[0]/10,cparray.shape[0],1):
        cparray[i,:,:] *= 0.0
    maxarray = numpy.amax(cparray, axis = 0)
    maxvalue = numpy.amax(maxarray[numpy.where(labelgrid == 1)])
    minvalue = numpy.amin(maxarray[numpy.where(labelgrid == 1)])
    return maxarray, minvalue, maxvalue


# catch indeces of the best curve
def best_curve_indeces(temp_maxarray, em_maxarray, labelgrid):
    temp_score = numpy.zeros(temp_maxarray.shape, dtype = int)
    em_score = numpy.zeros(em_maxarray.shape, dtype = int)
    for i in range (0,temp_maxarray.shape[0],1):
        for j in range (0,temp_maxarray.shape[1],1):
            val = temp_maxarray[i,j]
            temp_score[i,j] = abs(len(numpy.where(temp_maxarray*labelgrid > val)[0]) - len(numpy.where(temp_maxarray*labelgrid < val)[0]) + len(numpy.where(labelgrid == 0)[0]))
            val = em_maxarray[i,j]
            em_score[i,j] = abs(len(numpy.where(em_maxarray*labelgrid > val)[0]) - len(numpy.where(em_maxarray*labelgrid < val)[0]) + len(numpy.where(labelgrid == 0)[0]))
    return numpy.where((temp_score+em_score) == numpy.amin(temp_score+em_score))[0][0], numpy.where((temp_score+em_score) == numpy.amin(temp_score+em_score))[1][0]


# check if the URL source exists
def check_url(url):
    ret = requests.get(url)
    if ret.status_code == 200:
        return True
    else:
        return False


# check if the file exists in the URL, and figure out the filename and satellite's version
def find_goesfile(time):
    year = time[0:4]
    if ((time[0:4] in ['1996','1997','1998']) or ((time[0:4] == '1999') and (time[5:7] == '01') and (float(time[8:10]) < 15))):
        file_postfix = time[2:4]+time[5:7]+time[8:10]
    else:
        file_postfix = time[0:4]+time[5:7]+time[8:10]
    # checking the version of the satellite beginning from goes08
    gver = '00'
    if check_url('http://umbra.nascom.nasa.gov/goes/fits/'+year+'/go08'+file_postfix+'.fits') == True : gver = '08'
    if check_url('http://umbra.nascom.nasa.gov/goes/fits/'+year+'/go09'+file_postfix+'.fits') == True : gver = '09'
    if check_url('http://umbra.nascom.nasa.gov/goes/fits/'+year+'/go10'+file_postfix+'.fits') == True : gver = '10'
    if check_url('http://umbra.nascom.nasa.gov/goes/fits/'+year+'/go11'+file_postfix+'.fits') == True : gver = '11'
    if check_url('http://umbra.nascom.nasa.gov/goes/fits/'+year+'/go12'+file_postfix+'.fits') == True : gver = '12'
    if check_url('http://umbra.nascom.nasa.gov/goes/fits/'+year+'/go13'+file_postfix+'.fits') == True : gver = '13'
    if check_url('http://umbra.nascom.nasa.gov/goes/fits/'+year+'/go14'+file_postfix+'.fits') == True : gver = '14'
    if check_url('http://umbra.nascom.nasa.gov/goes/fits/'+year+'/go15'+file_postfix+'.fits') == True : gver = '15'
    filename = 'http://umbra.nascom.nasa.gov/goes/fits/'+year+'/go'+gver+file_postfix+'.fits'
    filename_short = './goes_rawdata/go'+gver+file_postfix+'.fits'
    return int(gver), filename, filename_short


# time conversion back to the datetime format
def converttime_sec_string(timing,start_time):
    stime = datetime.datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%S")
    stime = stime.replace(hour = 0)
    stime = stime.replace(minute = 0)
    stime = stime.replace(second = 0)
    if (isinstance(timing,float)):
        ctime = stime + datetime.timedelta(seconds = timing)
        stiming = datetime.datetime.strftime(ctime, "%Y-%m-%d %H:%M:%S")
    else:
        stiming = []
        for i in range (0, len(timing),1):
            ctime = stime + datetime.timedelta(seconds = timing[i])
            sctime = datetime.datetime.strftime(ctime, "%Y-%m-%d %H:%M:%S")
            stiming.append(sctime)
    return stiming


# correction for the negative flux values (instrumental)
def correct_flux(fluxes):
    # forward correction
    for i in range(1,fluxes.shape[0],1):
        if (fluxes[i,0] <= 0.0): fluxes[i,0] = fluxes[i-1,0]
        if (fluxes[i,1] <= 0.0): fluxes[i,1] = fluxes[i-1,1]
    # backward correction
    for i in range(fluxes.shape[0]-2,-1,-1):
        if (fluxes[i,0] <= 0.0): fluxes[i,0] = fluxes[i+1,0]
        if (fluxes[i,1] <= 0.0): fluxes[i,1] = fluxes[i+1,1]
    return fluxes

# apply Savitzky-Golay smoothing if selected
def sgsmooth_flux(fluxes, gver):
    if (gver <= '12'):
        window_size = 11
    else:
        window_size = 15
    order = 3
    fluxes[:,0] = savgol_filter(numpy.copy(fluxes[:,0]), window_size, order)
    fluxes[:,1] = savgol_filter(numpy.copy(fluxes[:,1]), window_size, order)
    return fluxes

# main
def TEBBS_calculate(start_time, end_time, plot_key = 0, sys_win = 0, savitzky_golay = 0, extend_50 = 0):

    # checking if the flare crossed the midnight point
    mid_cross = False
    if (start_time[0:10] != end_time[0:10]): mid_cross = True

    if (mid_cross == False):
        gver, goesfile, goesfile_short = find_goesfile(start_time)
        if (gver == 0): sys.exit("The corresponding GOES file was not found")
        if (os.path.isfile(goesfile_short) == False):
            ftemp = wget.download(goesfile)
            shutil.move(ftemp, './goes_rawdata/')
        timing, fluxes = read_flux(goesfile_short)
        flare_start_time = return_sec(start_time)
        flare_end_time = return_sec(end_time)

    if (mid_cross == True):
        gver1, goesfile1, goesfile_short1 = find_goesfile(start_time)
        gver2, goesfile2, goesfile_short2 = find_goesfile(end_time)
        if (gver1 == 0): sys.exit("The corresponding GOES file was not found")
        if (gver2 == 0): sys.exit("The corresponding GOES file was not found")
        if (gver1 != gver2): sys.exit("The versions of the files were inconsistent for the day transition")
        gver = gver1    # otherwise we've exit from the program
        if (os.path.isfile(goesfile_short1) == False):
            ftemp = wget.download(goesfile1)
            shutil.move(ftemp, './goes_rawdata/')
        if (os.path.isfile(goesfile_short2) == False):
            ftemp = wget.download(goesfile2)
            shutil.move(ftemp, './goes_rawdata/')
        timing1, fluxes1 = read_flux(goesfile_short1)
        timing2, fluxes2 = read_flux(goesfile_short2)
        timing2 += 86400.0
        timing = numpy.concatenate((timing1, timing2), axis = 0)
        fluxes = numpy.concatenate((fluxes1, fluxes2), axis = 0)
        # for windows-based distributives:
        if (sys_win == 1): fluxes = numpy.concatenate((fluxes2, fluxes1), axis = 0)
        flare_start_time = return_sec(start_time)
        flare_end_time = return_sec(end_time)+86400

    fluxes = correct_flux(fluxes)
    if (savitzky_golay == 1): fluxes = sgsmooth_flux(fluxes, gver)
    flare_peak_time = find_max_sec(timing, fluxes, flare_start_time, flare_end_time)
    if ((int(flare_peak_time) == 0) or (int(flare_peak_time) == int(flare_start_time))):
        print "The timing is incorrect. Interrupting run for the current flare..."
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0
    Amin, Bmin = extract_min(timing, flare_start_time, flare_peak_time, fluxes)
    bgrid = make_grid(Amin, Bmin)
    labelgrid = numpy.zeros([bgrid.shape[0],bgrid.shape[1]], dtype = int)+1
    fluxes_bsub, fluxes_bsub_ext, ftimimng, ftiming_ext = extract_fluxes(timing, flare_start_time, flare_peak_time, fluxes, bgrid)
    rising_phase_bins = fluxes_bsub.shape[0]

    # Now, one has a variety of background-subtracted arrays for both with first 1/6 of the rising phase
    # (called *ext) and without it. Let's calculate Temperatures and EMs for them.
    temp, em = calculate_temperature_em(fluxes_bsub, gver)
    temp_ext, em_ext = calculate_temperature_em(fluxes_bsub_ext, gver)
    # Now, let's apply a couple of tests. Linear fit suggested in the paper is not used for the tests.
    labelgrid, temperature_min_flag = test_temperature_min(temp, labelgrid)
    labelgrid = test_array_increase(temp, labelgrid)
    labelgrid = test_array_increase(em, labelgrid)
    labelgrid, initials_flag = test_t_em_initials(temp, temp_ext, em, em_ext, labelgrid)


    # Now, after all the tests performed, let's calculate T and EM and plot the variety of curves for them
    if (extend_50 == 0):
        flare_end_time += (flare_end_time-flare_start_time)*0.2
    else:
        flare_end_time += (flare_end_time-flare_start_time)*0.5
    flareflux, flareflux_ext, ftimimng, ftiming_ext = extract_fluxes(timing, flare_start_time, flare_end_time, fluxes, bgrid)
    plot_temp, plot_em = calculate_temperature_em(flareflux_ext, gver)
    temp_maxarray, temp_errmin, temp_errmax = error_range(plot_temp,labelgrid)
    em_maxarray, em_errmin, em_errmax = error_range(plot_em,labelgrid)

    if (sum(sum(labelgrid)) == 1):
        ibest,jbest = numpy.where(labelgrid == 1)
    else:
        ibest,jbest = best_curve_indeces(temp_maxarray, em_maxarray, labelgrid)

    # plotting the graph if the plot_key is set
    if (plot_key == 1):
        for i in range (0,temp.shape[1]):
            for j in range (0,temp.shape[2]):
                if (labelgrid[i,j] == 1):
                    plt.plot(ftiming_ext, plot_temp[:,i,j])
        plt.plot(ftiming_ext, plot_temp[:,ibest,jbest], linewidth = 5.0, color='grey')
        plt.show()
        for i in range (0,temp.shape[1]):
            for j in range (0,temp.shape[2]):
                if (labelgrid[i,j] == 1):
                    plt.plot(ftiming_ext, plot_em[:,i,j])
        plt.plot(ftiming_ext, plot_em[:,ibest,jbest], linewidth = 5.0, color='grey')
        plt.show()
        plt.imshow(labelgrid)
        plt.show()


    # calculating parameters
    if (extend_50 == 0):
        Tmax = numpy.amax(plot_temp[5*plot_temp.shape[0]/36:5*plot_temp.shape[0]/6,ibest,jbest])
        EMmax = numpy.amax(plot_em[5*plot_em.shape[0]/36:5*plot_em.shape[0]/6,ibest,jbest])
        Tmax_time = ftiming_ext[5*plot_temp.shape[0]/36 + numpy.argmax(plot_temp[5*plot_temp.shape[0]/36:5*plot_temp.shape[0]/6,ibest,jbest])]
        EMmax_time = ftiming_ext[5*plot_temp.shape[0]/36 + numpy.argmax(plot_em[5*plot_em.shape[0]/36:5*plot_em.shape[0]/6,ibest,jbest])]
    else:
        Tmax = numpy.amax(plot_temp[plot_temp.shape[0]/9:2*plot_temp.shape[0]/3,ibest,jbest])
        EMmax = numpy.amax(plot_em[plot_em.shape[0]/9:2*plot_em.shape[0]/3,ibest,jbest])
        Tmax_time = ftiming_ext[plot_temp.shape[0]/9 + numpy.argmax(plot_temp[plot_temp.shape[0]/9:2*plot_temp.shape[0]/3,ibest,jbest])]
        EMmax_time = ftiming_ext[plot_temp.shape[0]/9 + numpy.argmax(plot_em[plot_em.shape[0]/9:2*plot_em.shape[0]/3,ibest,jbest])]
    # converting the times for flares
    ftiming_ext = converttime_sec_string(ftiming_ext,start_time)
    Tmax_time = converttime_sec_string(Tmax_time, start_time)
    EMmax_time = converttime_sec_string(EMmax_time, start_time)
    # creating output arrays of T and EM
    T_out = plot_temp[:,ibest,jbest]
    EM_out = plot_em[:,ibest,jbest]
    fluxes_out = flareflux_ext[:,ibest,jbest,:]
    # Figuring out the peaks and the peak times of the fluxes
    if (extend_50 == 0):
        flxlength = len(numpy.copy(fluxes_out[:,0]))
        AFluxMax = numpy.amax(numpy.copy(fluxes_out[5*flxlength/36:5*flxlength/6,0]))
        BFluxMax = numpy.amax(numpy.copy(fluxes_out[5*flxlength/36:5*flxlength/6,1]))
        AFluxTime = ftiming_ext[5*flxlength/36 + numpy.argmax(numpy.copy(fluxes_out[5*flxlength/36:5*flxlength/6,0]))]
        BFluxTime = ftiming_ext[5*flxlength/36 + numpy.argmax(numpy.copy(fluxes_out[5*flxlength/36:5*flxlength/6,1]))]
    else:
        if (len(fluxes_out.shape) > 2):
            print "ALERT: array reformation. Check the results."
            fluxes_out = numpy.copy(fluxes_out[:,0,:])
        flxlength = len(numpy.copy(fluxes_out[:,0]))
        AFluxMax = numpy.amax(numpy.copy(fluxes_out[flxlength/9:2*flxlength/3,0]))
        BFluxMax = numpy.amax(numpy.copy(fluxes_out[flxlength/9:2*flxlength/3,1]))
        AFluxTime = ftiming_ext[flxlength/9 + numpy.argmax(numpy.copy(fluxes_out[flxlength/9:2*flxlength/3,0]))]
        BFluxTime = ftiming_ext[flxlength/9 + numpy.argmax(numpy.copy(fluxes_out[flxlength/9:2*flxlength/3,1]))]
    return fluxes_out, T_out, EM_out, ftiming_ext, AFluxMax, AFluxTime, BFluxMax, BFluxTime, Tmax, temp_errmin, temp_errmax, Tmax_time, EMmax, em_errmin, em_errmax, EMmax_time, temperature_min_flag, initials_flag, rising_phase_bins
