![http://ceesdekkerlab.tudelft.nl/wp-content/uploads/GUI_detect.png](http://ceesdekkerlab.tudelft.nl/wp-content/uploads/GUI_detect.png)

# Quick manual version #

  1. Select file type, filter frequency, and voltage.
  1. Open the directory.
  1. Double click on a trace file and/or select a segment.
  1. Click **Find Events** and check that events found are not junk.
  1. Click **Add to Master Trace**.
  1. Repeat steps 3 to 5 for all of the traces you want to analyze.
  1. Click **Save** and the program will export all detected events.

# Quick automated version #

  1. Sort all of your experiments into separate folders.
  1. Make sure any non-trace files or folders have a name starting with a tilde ~ so that the script will ignore them.
  1. Open each folder and select your detection parameters for that particular folder. For each folder check your detection parameters by finding events in several traces. Once you are happy click **Create PAR**.
  1. Repeat step 2 for all of the different folders.
  1. Go to the root directory containing all of these folders and check **Auto Save** and **Cycle Subfolders**.
  1. Click **Cycle All**.

# Detailed guide #

Now the files are ready to analyze. Open the Detect Events script (GUI\_detect) and open the folder containing the traces you want to analyze.
  * Select the proper file type
    * LabView Binary
    * Clampex Text (.atf)
    * Text File
    * Current Only Text File (**Note:** here you must manually set the timestep)
    * ABF Axon Binary (.abf)
  * In most experiments the SNR is low enough that the traces must be digitally filtered. Turn on filtering and select Gaussian (minimal distortion). Enter your filtering frequency.
  * The box immediately to the right of the filtering frequency is the gain a\*B. This is only used with certain file types to have the proper current levels. This is not used for ABF Axon Binary files.
  * The voltage box (mV) is used to calculate the conductance values in nS, the proper measurement voltage for the files in the current folder should be set.
  * The resistance box below the voltage shows the resistance calculated using the resistance(Mohms) = voltage(mV)/current(nA). The resistance is recalculated everytime you change the voltage or load a new trace.
  * The save unfiltered check box creates a second trace file for the events with the unfiltered data. This is mostly unused except in the Ghost feature of GUI\_events. If it is checked the results will take longer to save. **Note:** Because of delay introduced by the filtering, the filtered and unfiltered traces are slightly shifted with respect to each other.
  * Set your detection parameters:
    * MovAvg Window: this is the size of the moving average window in points. It should be smaller for unstable baselines and larger for stable baselines. Ideally the baseline should not change much within the size of this window.
    * Peak Det: The peak detection factor (PDF). This sets the detection limit using = PDF\*STD. This is how far away from the moving average line an event must be in order for it not to be treated as noise. This level is shown as the red line on the trace plot. Almost always, this is kept at a value of 5. If you go to low it will start to pickup noise. If you go to high you may miss translocation events.
    * Min Dwell: The minimum dwell time (ms) for an event to be accepted.
    * Max Dwell: The maximum dwell time (ms) for an event to be accepted.
    * MinEventsGap: The minimum number of points in between two successive events. This should be left as two unless you know what you're doing.
    * Extra Points: The number of extra points on either side of the event to save when extracting the events from the trace files. Typically a value of 100 is used. (This will also give you your baseline peak in GUI\_events’s current histogram).
    * Event Type: At the bottom of the detection parameter box. Choose from:
      * Down events (normal setting, high salt)
      * Up events (low salt experiments)
      * Both (intermediate salt levels [Kowalczyk SW, Dekker C. Measurement of the docking time of a DNA molecule onto a solid-state nanopore](see.md) **Note:** This option uses an algorithm which is much slower, than the other two options above and has not been thoroughly tested yet!
    * Iterative Detection: This applies the same detection multiple times to the same file and each time it replaces the duration of detected events in the trace file with the value of the baseline at the start of the event. Then it recalculates the moving average. Because events have effectively been removed, distortions in the moving average are also removed. This function can be useful if you have events very close together or the event rate is high. The more iterations you run the more it should converge to the proper number of detected events (if all parameters are set correctly).

  * Baseline (how should we determine the local baseline):
    * Moving Average: The standard option used in analysis.
    * Manual (Initial Guess): Enter a value for the baseline in the text box. The first iteration of the detection will use this value. Subsequent loops will calculate the moving average (on the modified trace).
    * Manual (Locked): Always use the value in the text box as the value of the baseline.
    * Current Hist Max: Take the maximum value of the current histogram as the baseline value. Set the number of bins such that the resolution is sufficient to provide a good estimate. **Note:** This will not work if there are more points within the events than the baseline in a given trace. The current histogram button will open a plot showing the peaks. This should be used to set the proper binning. ![http://ceesdekkerlab.tudelft.nl/wp-content/uploads/TR_current.png](http://ceesdekkerlab.tudelft.nl/wp-content/uploads/TR_current.png)
  * Sigma (How to detect the value of the Standard Deviation / Noise Level):
    * Trace STD: Take the standard deviation of the entire trace. This only works if you have a very stable baseline.
    * Manual (Initial Guess): Enter a value for sigma in the text box. The first loop of the detection will use this value. Subsequent loops will calculate sigma on the global trace (on the modified trace) **Note:** This is to be used for a trace with many events but stable baseline.
    * Manual (Locked): Always use the value in the text box as the value of sigma.
    * Moving STD Hist Peak: Calculate STD over a local window with size in points given in the text box below. Take the maximum point of the local STD histogram. The STD Hist button will open a plot with these peaks. This should be used to set the proper binning. ![http://ceesdekkerlab.tudelft.nl/wp-content/uploads/TR_std.png](http://ceesdekkerlab.tudelft.nl/wp-content/uploads/TR_std.png)
    * Moving STD: NOT IMPLEMENTED YET!

  * Once your detection parameters are set you should save them by clicking **Create PAR**. This will create a file called ~detect.par which contains all of the detection settings you just chose. To reload these detection settings at a later date, use the **Load PAR** button.

  * Next you can double click a trace file to plot it. The magenta line at the center of the trace will show the moving average, with the two green lines showing the moving average +/- sigma. The red line shows the detection level. These lines can be removed by unselecting **Diagnostic Plot**. At any time you can open this plot externally by selecting **Plot Trace Win**.
  * Clicking **Find Events** will then detect all events using the specified parameters.
  * Once events are found you can browse through them with the Event Info panel buttons **Up** and **Down**. A variety of event information is displayed. The plot of the event can be opened externally by clicking **Plot Event**.
  * In the trace plot events will highlighted in red.
  * Once you are satisfied with the detected events, you can click **Add to Master Trace**. This creates a trace file with all of the detected events concatenated together. At any time you can browse either the events detected in the current trace file loaded (local) or the events detected since the start of the analysis (master). To plot the master trace file click either **Plot Master Win** (external plot) or **Plot Master Here** (in GUI).
  * This process can be repeated until you have sufficient events. The events in the master trace can be saved at any time by clicking **Save**. This creates a new directory called ~transalyzed-10000Hz. Where a filtering frequency of 10000Hz was used in this example.
  * To export the events to OpenNanopore select OpenNanopore files. This will create a directory called OpenNanopore within the ~transalyzed-10000Hz directory once you save the master trace. There will be two files:
    * **events.mat** – a trace file containing all of the events
    * **OpenNanopore\_parameters.txt** – the parameters to be used in the detection
  * You can analyze all of the trace files in a directory by clicking **Cycle All**. If the **Auto Save** option is turned on it will automatically save the master file at the end, otherwise you will have to do it manually by clicking Save. **Note:** **Cycle All** is the same as manually loading each file, clicking **Find Events** and then **Add to Master Trace**.

Analysis can be automated across multiple experiments (each in its own folder) by first creating PAR files in each folder to be analyzed, then going to the root folder and selecting **Cycle All** while the **Auto Save** and **Cycle Subfolders** options are on. **Note:** Folders not containing a PAR file will not be analyzed.

**Note:** Clicking **Dump Variables** will save all of the working variables to the current Matlab workspace.

These are the files created when you click **Save**:

**~detect.par**
A copy of the PAR file used.

**analysis\_events**
Each detected event is a new row with the following columns:
```
Column 1 Avarage Amplitude with local minima method (nA)
Column 2 Maximum Amplitude for the event (nA)
Column 3 Integral of event in nA*ms (pC)
Column 4 FWHM Dwell time (ms)
Column 5 Event Baseline (nA)
Column 6 Detection Level (nA)
Column 7 Event Type (not used?)
Column 8 start point index in rawtrace 
Column 9 end point index in rawtrace  
```

**analysis\_extra\_info.txt**
```
Column 1 file number 
Column 2 Resistance for this trace file (Mohms)
Column 3 Total time of this local trace (s)
Column 4 Total time passed in entire experiment (s)
Column 5 Number of local minima in event 
Column 6 Location of first local minima
Column 7 Location of last local minima
Column 8 Dwell time calculated using xend – xstart (ms)
Column 9 Amplitude in nA, calculated by Integral/FWHM_Dwell_time
```

**analysis\_events\_min\_max**
```
Column 1 min current level inside event
Column 2 max current level inside event
Column 3 min current level in extrapoints
Column 4 max current level in extrapoints
Column 5 std of current in extrapoints
Column 6 mean current level in extrapoints
Column 7 std of current in extrapoints at start (before translocation)
Column 8 mean current level in extrapoints at start (before translocation)
Column 9 std of current in extrapoints at end (after translocation)
Column 10 mean current level in extrapoints at end (after translocation)
```

**analysis\_files\_info.txt**
Column 1 contains the filename in quotes, Column 2 contains the filenumber (same as used in (analysis\_extra\_info.txt)


**analysis\_parameters.txt**
The various parameters used in the analysis, legacy file, use ~detect.par instead.

**analysis\_rawtrace**
The Master Trace, one column with all of the filtered current points in nA.

**analysis\_unfiltered\_rawtrace**
The Master Trace, one column with all of the current points in nA. No filter applied.

**analysis\_rawtrace\_baseline**
The Master Trace, one column with all of the filtered current points in nA. Here all event’s baselines are shifted to zero.