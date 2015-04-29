![http://ceesdekkerlab.tudelft.nl/wp-content/uploads/GUI_events.png](http://ceesdekkerlab.tudelft.nl/wp-content/uploads/GUI_events.png)

Now we can go and do statistics on the detected events. Open the GUI\_events m-file and run it.
**Note:** By default the script will open the rawtrace file. If your rawtrace file is too large Matlab may run out of memory! If you don't want to load the rawtrace unselect Load Rawtrace too in the top right corner.

In the top right there are two drop down selection boxes which define which parameters to use for the analysis. These should be selected before you load the data:
  * Dwell Time (how to calculate the events dwell time):
    * New Dwell (FWHM) – this option should be always used. It corresponds to the Full Width Half Maximum dwell time.
    * Old Dwell (t\_end – t\_start) – an old legacy variable which defines the dwell time using the points where the trace crosses the baseline at the start and end of the event. This is highly inaccurate and very prone to error.
  * Amplitude:
    * Avg Amplitude (bet. loc. Minima): The average level between the first and last local minima in the event. This is prone to error due to long tails having a local minimum. Do not use.
    * Max Amplitude: The maximum point in the event. To be used for proteins, short DNA, and other short spikes. The value here can be highly dependent on the filtering frequency especially at short dwell times.
    * Integral/FWHM: This has been found to be the most accurate but has the downside that it can produce unclear population distributions if you have various folded DNA events. Use caution.

Click **Load Events File** and select the **analysis\_events** file you want to open. The files will be loaded into Matlab and if all goes well you will see a scatter plot in the plot on the bottom right and a list of the events detected in the main listbox.

If you have loaded the rawtrace file you can click on an event in this list and it will be plotted in the graph on the bottom left. Selecting **Show Baseline** will show the baseline level as a green line in this plot. If you saved the unfiltered data in GUI\_detect you can overlay the unfiltered trace by selecting **Overlay Unfiltered**. You can save this plot by clicking **Open Event Plot**. You can also save this event's trace points by clicking **Save Event Dat**. This will create a file in /plots called event\_2\_filtered.dat where 2 is the number of the event in this example. The Ghost function will be described in a later version.

The voltage (mV) level is automatically loaded. Make sure it is correct if you want to use units of conductance (nS) for your analysis.

Events can be selected through a large number of parameters as listed in the Event Selection (Salt) Criteria box on the right. Possible min-max selection parameters are:
  * Dwell time in ms
  * Amplitude. **Note:** this selection is always done in nA regardless of whether you have the Conductance option selected.
  * Baseline: **Note:** this selection is always done in nA regardless of whether you have the Conductance option selected. **Note:** This is useful for removing clogs.
  * Event Number: Useful if end of experiment was bad or if there are too many events.
  * ECD Integral: Useful to sort by integrated area.
  * Use Fold Counts: to be described later
  * FEXPtMnMax and LEXPtMnMax: The maximum allowed deviation of the mean of the first and last extra points relative to a 0 baseline level.

To preview which population of events you will select you can click the **Draw** button and this will show the chosen limits on the plot. Once you have chosen the proper selection values click **Select Events**. To save these selection parameters and use them at a later date there are two options:
  * The Save Cln button will create a file called clean.salt with the selection parameters which can be reloaded later using the **Load Cln** button. This is typically used to select a population of events or remove the outliers and clogs from the dataset.
  * If you want to choose your own filename use the **Save Salt** / **Load Salt** buttons.

There are nine possible plot types as shown in the Plot 2 Type box on the bottom right corner:
  * Scatter: Plots Amplitude vs. Dwell Time. This plot can have a log scale for the Dwell times by selecting Dwell X Scale Log.
  * Baseline: This plots the value of the baseline at each detected event. Any clogs present in the data should be evident as large drops in the value of the baseline.
  * Rpore: The resistance of the pore as a function of experimental time.
  * Dwell: The dwell time histogram. The number of bins can be chosen with the text box on the right.
  * CB: Conductance blockade histogram. The number of bins can be chosen with the text box on the right.
  * ScatterInt: The integral of each event versus its dwell time.
  * Current: A histogram of all the current points in each event. The baseline for each event is shifted to zero. **Note:** This only works when the rawtrace is also loaded. Also if you select **Use Extra Points** it will include the extra points on either side of each event and your baseline peak will be better defined.
  * IntHist: The histogram of all the integral values. The number of bins can be chosen with the text box on the right.
  * HeatMap: A scatter plot showing the density in each area. This is useful for large datasets.

For each of these plots you can click Plot to open it in an external window.

Clicking **Mass Plot Export** will save all of the plots to the \plots directory in fig, eps, large png, and small png formats as well as the raw data used in each plot.

**Note:** Clicking **Dump Variables** will save all of the working variables to the current Matlab workspace.