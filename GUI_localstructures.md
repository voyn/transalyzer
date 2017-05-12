![GUI_localstructures](../../raw/wiki/GUI_localstructures.png)

This GUI can be used to find and characterize spikes present within translocation events.

The following describes the typical procedure for detection of spikes within events. Note that peaks within events can be due to noise, folding, DNA knots, co-translocation of multiple molecules.

* In GUI_events select a population of DNA events you wish to analyse. In general, for the analysis. Try to remve super short events.
* In GUI_events use the current histogram to determine the level of the dsDNA blockade (I_1) for your experiment. Note: the program assumes no significant changes in the value of I_1 over the course of the experiment.
* Start GUI_localstructures and enter the value of the DNA blockade (I_1) into the Baseline text box.
* Choose "Open Events" and select an analysis_events file to analyse.
* Click "Region Analysis" in the "1. Calc selection param" box. This calculates the Normalized Integral between I_1 and I_2 for all events, as described in the publication.
* A histogram of the Normalized Integral between I_1 and I_2 will be plotted. Values close to 1 represent fully folded events, while values close to 0 are unfolded events. Those in-between represent partially folded event or events with spikes. Typically short duration spikes do not contribute much to the normalized integral.
* In the "2. Classify Events" box, select the highest Normalized Integral event to be analysed under "Region Int". For events with short spikes we often set the Max limit between 0.2 and 0.4. The goal here is the filter out highly folded events from the data set.
* Click select.
* Click "Assign" to classify these events as "Unfolded - Unsorted".
* In box "4. Peak Finder" select a Threshold and Selection factor. The Threshold sets the minimum amplitude for a peak to be detected in units of I_1 relative to the open pore current I_0. So a Threshold of 2.5 will mean peak's over 2.5xI_1 from the open pore baseline will be detected. The Selection factor sets how far above the surrounding data a peak must be to be counted. See documentaion for the peakfinder function: [https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder-x0--sel--thresh--extrema--includeendpoints--interpolate-](https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder-x0--sel--thresh--extrema--includeendpoints--interpolate-)
* Click "Detect Local Structures". This detects all peaks using the given settings in all events.
* In box "5. Peak Stats and Plot Gen Panel", under Analyse Peaks from select "Unfolded - All". This will plot only a subset of the peaks based on how events have been classified.
* Now a number of plots can be generated:
  * Scatter Spikes creates a scatter plot of peak amplitude versus peak dwell time.
  * Dwell hist - is a histogram of peak durations.
  * Amp hist - is a histogram of peak amplitudes.
  * Pos (Norm) hist - shows a histogram of the normalized center position of the peaks.
  * Pos (Abs) hist - shows the absolute position of the detected peaks.
  * Current hist - is the same as in GUI_events but using only the events which have been selected.
  * Peaks/Events - shows a histogram of the number of peaks per event.
* "Mass Plot" will export all plots at once.
* Time Hist - will show the fraction of events with peaks as a function of experimental time.
* The "Dump Variables" button will place all analysis variables into the Matlab base workspace. All variables are the same as in GUI_events except for:
  * event_local_struct_info
  * here each event is a row and the columns are:
    * 1 event id - a unique identifying number for each event
    * 2 integral - the event charge deficit (ECD)
    * 3 FWHM_dwell - duration of the event (ms)
    * 4 number of local struct - the number of peaks detected
    * 5 event start ind - index of event start in rawtrace_vect
    * 6 event end ind - index of event end in rawtrace_vect
    * 7 region integral - the normalized integral between I_1 and I_2
    * 8 max peak FWHM - duration of the largest peak
    * 9 group 1 event classification type:
      * value -1 - bad
      * value 0 - unsorted
      * value 1 - folded
      * value 2 - clogged
      * value 3 - unfolded
    * 10 group 2 event classification type:
      * value 0 - unsorted
      * value 1 - bare (no spikes)
      * value 2 - spike A
      * value 3 - spike B
      * value 4 - spike C
      * value 5 - spike D
      * value 6 - spike E
      * value 7 - spike F
      * value 8 - spike G
      * value 9 - spike H
      * value 10 - spike I
    * 11 number of local struct from the region analysis
    * 12 mean value extra points before event
    * 13 mean value extra points after event
    * 14 event selected: 1 -yes 0 - no
  * peak_info
    * each row is a peak and the columns are:
    * 1 - eventid - ID of the event with this peak
    * 2 - peakid - unique ID for this peak
    * 3 - peak position index - index of the peak maximum
    * 4 - peak position relative to FWHM start in time - time between event start and peak (ms)
    * 5 - peak FWHM - peak duration
    * 6 - peak FWHM start time
    * 7 - peak FWHM start mag
    * 8 - peak amplitude in nA not including first DNA baseline (I_1)
    * 9 - peak FWHM end magnitude
    * 10 - peak position normalized with event FWHM

