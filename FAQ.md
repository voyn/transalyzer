#Frequently Asked Questions

# Frequently Asked Questions #

## 1. How can I open my file type with Transalyzer? ##

If your data is stored in a file type other than those currently supported, you must modify the fun\_load\_trace\_file function in GUI\_detect. Add a new case at the end of the switch statement:

```
        % get the trace file type
        file_type = get(handles.pop_file_type,'Value');
                
        switch  file_type % switch based on type of trace
                    
             case # % add a new case after last case statement
```

Your code may use two inputs:
  1. filename - complete path to trace file
  1. currentsegment - which segment of the trace to open

And generate the following vectors and variables:
  1. handles.file\_load\_successful - file opened correctly 1, or not 0
  1. handles.trace - current values of the trace in nA
  1. handles.time\_vector - timepoints for each current point
  1. handles.timestep - time between each current point
  1. handles.total\_number\_of\_segments - total number of segments in this trace

Finally, you must open GUI\_detect.fig using [guide](http://www.mathworks.nl/help/matlab/ref/guide.html) and add your file type to the list of types in the drop down menu (_handles.pop\_file\_type_).

## 2. How do I diagnose a programming issue? ##

In Matlab GUIs, data is typically stored in structures called [handles](http://www.mathworks.nl/help/matlab/ref/guidata.html), which can not be directly accessed from the base workspace. One useful method for debugging is to pass the variables used into the base workspace just before the point where the error occurs. This can be done using the [assignin](http://www.mathworks.nl/help/matlab/ref/assignin.html) function:
```
assignin('base', 'name_of_variable_in_base_workspace', handles.name_of_variable_in_GUI)
```