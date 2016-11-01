from TEBBS import *

# example of usage of TEBBS and printing out peak T with its error bounds
start_time = '2014-06-13T06:43:00'
end_time = '2014-06-13T06:51:00'
fluxes, Tarray, EMarray, timing, T, Tmin, Tmax, Ttime, EM, EMmin, EMmax, EMtime = TEBBS_calculate(start_time, end_time, plot_key = 0)
print T, Tmin, Tmax