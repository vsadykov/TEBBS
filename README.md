# TEBBS
The algorithm proposed in the paper by Ryan et al. (2012) "The Thermal Properties of Solar Flares over Three Solar Cycles Using GOES X-Ray Observations"

### Input parameters:

start_time - start time of the flare in 'yyyy-mm-ddThh:mm:ss' format

end_time - end time of the flare in 'yyyy-mm-ddThh:mm:ss' format

### Optional:

plot_key - key to show the plot via local pyplot interface. 0 by default ( = do not show)

### Output parameters:

fluxes - fluxes of the GOES X-ray data during the flare

Tarray - best solution for the Temperature curve

EMarray - best solution for the Emission Measure curve

timing - timing of the flux in 'yyyy-mm-dd hh:mm:ss' format

T, Tmin, Tmax - peak value of the Temperature and its error bounds

EM, EMmin, EMmax - peak values of the Emission Measure and its error bounds

tmin_test_flag - flag showing if the hot flare test passed (1 for passed, 0 for failed).
If the test is failed, it is not applied for the studied event.

init_test_flag - flag showing if the test for the initial 1/6 interval of the growing phase passed (1 for passed, 0 for failed).
If the test is failed, the curve having the smallest oscillation of initials is chosen.

### Example:

start_time = '2014-06-13T06:43:00'

end_time = '2014-06-13T06:51:00'

fluxes, Tarray, EMarray, timing, T, Tmin, Tmax, Ttime, EM, EMmin, EMmax, EMtime, tmin_test_flag, init_test_flag = TEBBS_calculate(start_time, end_time, plot_key = 0)
