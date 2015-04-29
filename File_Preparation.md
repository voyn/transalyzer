# Sort your data #

Move all IV curves and noise traces into a new folder. Move each of the traces corresponding to a different voltage or experiment into its own folder as shown below.

**Note**: The program will ignore any files or folders if their name starts with a tilde ~

![http://ceesdekkerlab.tudelft.nl/wp-content/uploads/TR_folder_screenshot.png](http://ceesdekkerlab.tudelft.nl/wp-content/uploads/TR_folder_screenshot.png)

If you are recording using LabView, set each DATALOG file to be around 10Mb at 500kS/s acquisition. Before starting analysis, the experimental data must be in the proper (temporal) order when you look at the directory in Matlab.

If your files end in numbers, make sure they use sufficient leading zeros to have the proper file order. If you need to add some leading zeros you can use the following procedure:

Install http://www.bulkrenameutility.co.uk/ while logged in as localadmin. Open the folder with your data (note: you can load data from subfolders under the Selections (12) box, check subfolders and files, uncheck folders).

If you have any periods in your filenames, these must be removed. The following procedure will convert the periods to underscores. Select all of your trace files. Check RegEx(1) and Extension(11). For Extension(11) select Remove. For RegEx(1) check Include Ext. and use
```
(.+)[.](.+)$
```
for Match and
```
\1_\2
```
for Replace. Run this once for every period in your filename. So if you have some file names containing three periods click Rename three times, each time the period closest to the end of the filename will be replaced.

![http://ceesdekkerlab.tudelft.nl/wp-content/uploads/TR_bulk1.png](http://ceesdekkerlab.tudelft.nl/wp-content/uploads/TR_bulk1.png)

Once your files no longer contain periods in their names, select all of your trace files and uncheck Extension(11), check RegEx(1) and uncheck the Include Ext option. The new file names will be shown on the right side of the file list. Files that will have their name changed will be colored green. Three RegEx(1) filters will need to be run:

For single digits (will change 1 into 0001)
```
Match:
(.+)\s(\d)$
Replace:
\1 000\2
```

For two digits (will change 12 into 0012)
```
Match:
(.+)\s(\d\d)$
Replace:
\1 00\2
```


For three digits (will change 123 into 0123)
```
Match:
(.+)\s(\d\d\d)$
Replace:
\1 0\2
```

![http://ceesdekkerlab.tudelft.nl/wp-content/uploads/TR_bulk2.png](http://ceesdekkerlab.tudelft.nl/wp-content/uploads/TR_bulk2.png)