![http://ceesdekkerlab.tudelft.nl/wp-content/uploads/TR_gaussian.png](http://ceesdekkerlab.tudelft.nl/wp-content/uploads/TR_gaussian.png)

If you are interested in determining the peaks (Gaussian fits) of dwell, current, or amplitude histograms you can click on Gaussian Fit to open up the GUI\_gaussian script.
  * First select the Data Input type:
    * Dwell
    * Blockade
    * Current
  * Next determine how you want to fit:
    * Single Gaussin (histfit): This will do a single fit.
    * Multi Gaussian (GMM script): This will try to fit multiple Gaussians to the data. The number of Gaussians to be fit is given in the # of Gaussians text box. The fit parameters of the first three distributions will be displayed.
    * Manual Gaussian (nofit): You can manually plot distributions, or more usefully you can adjust previously fitted distributions.

Once you are ready click **Replot**. The button **Open Plot** will open the plot in an external window where it can be saved.