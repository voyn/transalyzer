function varargout = GUI_localstructures(varargin)
% GUI_LOCALSTRUCTURES M-file for GUI_localstructures.fig
%      GUI_LOCALSTRUCTURES, by itself, creates a new GUI_LOCALSTRUCTURES or raises the existing
%      singleton*.
%
%      H = GUI_LOCALSTRUCTURES returns the handle to a new GUI_LOCALSTRUCTURES or the handle to
%      the existing singleton*.
%
%      GUI_LOCALSTRUCTURES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_LOCALSTRUCTURES.M with the given input arguments.
%
%      GUI_LOCALSTRUCTURES('Property','Value',...) creates a new GUI_LOCALSTRUCTURES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_localstructures_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_localstructures_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help GUI_localstructures

% Last Modified by GUIDE v2.5 29-Jul-2014 16:13:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_localstructures_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_localstructures_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GUI_localstructures_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_localstructures (see VARARGIN)

% Choose default command line output for GUI_localstructures
handles.output = hObject;

setappdata(0,'HandleMainGUI',hObject);

% Update handles structure
guidata(hObject, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = GUI_localstructures_OutputFcn(~, ~, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_exit_Callback(~, ~, ~)

    delete(gcf)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this code is from GUI_events button_load_events_file_Callback
function button_open_Callback(hObject, ~, handles)

    clear rawtrace_vect

    set(handles.edit_status,'String','Waiting for file...')
    [handles.filename, handles.pathname] = uigetfile( ...
    {  '*',  'All Files (*.*)'; ...
       '*_events','Events File'}, ...
       'Pick a file'); % get the events file

    set(handles.listbox_eventsbox,'Value',1); % in case old selection value higher than new total (which causes events box to disappear)

    fun_load_analysis_files(handles.pathname)
    
    handles = guidata(gcbo);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update various text boxes / radio buttons with the new data, and/or default values and
    % states
    set(handles.edit_status,'String','Opening events file...done')
    
    % local struct code
    event_id = 1:1:size(handles.event_info,1);
    event_local_struct_number = zeros(size(handles.event_info(:,2)));
    
    handles.event_local_struct_info = [event_id' handles.event_info(:,3) handles.event_info(:,4) event_local_struct_number event_local_struct_number event_local_struct_number event_local_struct_number event_local_struct_number event_local_struct_number event_local_struct_number event_local_struct_number event_local_struct_number];
    % 1 id
    % 2 integral
    % 3 FWHM_dwell
    % 4 number of local struct
    % 5 event start ind
    % 6 event end ind
    % 7 region integral
    % 8 max peak FWHM
    % 9 group 1 type:
    %       -1 - bad
    %       0 - unsorted
    %       1 - folded
    %       2 - clogged
    %       3 - unfolded
    %       
    % 10 group 2 type
    %       0 - unsorted
    %       1 - bare (no spikes)
    %       2 - spike A
    %       3 - spike B
    %       4 - spike C
    %       5 - spike D
    %       6 - spike E
    %       7 - spike F
    %       8 - spike G
    %       9 - spike H
    %       10 - spike I
    % 11 number of local struct from region analysis
    % 12 mean first extra points
    % 13 mean last extra points
    % 14 selected 1 -yes 0 - no
    
    handles.peak_info = [];
    % 1 - eventid
    % 2 - peakid
    % 3 - peakpos_ind
    % 4 - peakpos_rel to FWHM start in time
    % 5 - peak FWHM
    % 6 - peak FWHM start time
    % 7 - peak FWHM start mag
    % 8 - peak amp nA not including first DNA baseline
    % 9 - peak FWHM end mag
    % 10 - peakpos normalized with event FWHM
    
    handles.region_analysis_peak_info = [];
    % 1 - eventid
    % 2 - peakid
    % 3 - peakpos_ind
    % 4 - peakpos_rel to FWHM start in time
    % 5 - peak FWHM   
    % 6 - peak FWHM start time
    % 7 - peak FWHM start mag
    % 8 - peak amp nA   
    % 9 - peak FWHM end mag
    % 10 - peakpos normalized with event FWHM
    
    % handles.event_more_info_temp = handles.event_more_info; % select everything
    % this contains :
    % 1 file number
    % 2 resistance
    % 3 time - local
    % 4 time - global
    % 5 number of local min
    % 6 index first local min
    % 7 index last local min
    % 8 dwellold (between start and stop)
    % 9 amplitude good
    
    guidata(hObject, handles); % update handles
    
    fun_update_listbox()
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_load_analysis_files(location)

    handles = guidata(gcbo);
    handles.pathname = location;
    handles.filename = 'analysis_events';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.edit_status,'String','Opening events file...')
    fid = fopen(strcat(handles.pathname,filesep,handles.filename)); % open events file
    if fid == -1
        disp(strcat(handles.filename,' could not be loaded!!')) % not good
    else
        handles.event_info = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f')); % grab the event info
    %                 1 amplitude_minima_maxima
    %                 2 maximum	Amplitude
    %                 3 integral
    %                 4 FWHM_dwell
    %                 5 baseline
    %                 6 detectionlevel
    %                 7 eventtypeUpDown	down =0 up=1
    %                 8 startpoint
    %                 9 endpoint
        if get(handles.check_flip,'Value')
            handles.event_info(:,1) = -1.*handles.event_info(:,1);
            handles.event_info(:,2) = -1.*handles.event_info(:,2);
            handles.event_info(:,3) = -1.*handles.event_info(:,3);
            handles.event_info(:,5) = -1.*handles.event_info(:,5);
            handles.event_info(:,6) = -1.*handles.event_info(:,6);
        end
        fclose(fid);
    end
    
    set(handles.edit_status,'String','Opening events file...done...opening rawtrace file...')

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filenameparam = strcat(handles.pathname, filesep,'analysis_parameters.txt');
    if isequal(exist(filenameparam,'file'),2) % if the paramter file exists (legacy)
        fid = fopen(filenameparam);
        
        if fid == -1
            disp(strcat(filenameparam,' could not be loaded!!'))
        else % file loaded successfully
            line = fgetl(fid); % skip Transalyzer Analysis
            line = fgetl(fid); % skip version
            line = fgetl(fid); % skip Analysis Date
            line = fgetl(fid); % date, not used anywhere
            line = fgetl(fid); % skip
            handles.voltage = str2double(fgetl(fid)); % Voltage (mV)
            line = fgetl(fid); % skip
            handles.fsample = str2double(fgetl(fid)); % Sampling Frequency (Hz)
            handles.timestep = 1/handles.fsample;
            line = fgetl(fid); % skip
            handles.localtracetime = str2double(fgetl(fid)); % Local Trace Time (s)
            line = fgetl(fid); % skip
            handles.peak_detection_factor = str2double(fgetl(fid)); % Peak Det Factor
            line = fgetl(fid); % skip
            handles.sigma_from_detect_GUI = str2double(fgetl(fid)); % Sigma (nA)
            line = fgetl(fid); % skip
            handles.readextrapoints = str2double(fgetl(fid)); % Extra Points (points)
            line = fgetl(fid); % skip
            handles.totaleventsfound = str2double(fgetl(fid)); % Total Events Found
            line = fgetl(fid); % skip
            handles.totaltimetrace = str2double(fgetl(fid)); % Total Time
            
            fclose(fid); % close the file
        end
        
        
    end

    % handles.event_info contains everything

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [fileactual filediscard] = strread(handles.filename, '%s %s', 'delimiter','_'); % legacy code, before standardized naming
    handles.filenametrace = char(strcat(fileactual,'_rawtrace')); % determine rawtrace file name
    fid = fopen(strcat(handles.pathname,filesep,handles.filenametrace));
    if fid == -1
            disp(strcat(handles.filename,' could not be loaded!!'))
    else
        handles.rawtrace_vect = cell2mat(textscan(fid, '%f')); % read in rawtrace file
        if get(handles.check_flip,'Value')
            handles.rawtrace_vect = -1.*handles.rawtrace_vect;
        end
        fclose(fid); % close the file
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [fileactual filediscard] = strread(handles.filename, '%s %s', 'delimiter','_'); % legacy code, before standardized naming
    handles.filenametrace_unfilt = char(strcat(fileactual,'_unfiltered_rawtrace')); % determine rawtrace file name
    if isequal(exist(handles.filenametrace_unfilt),2) % if the extra info file exists (legacy)
        fid = fopen(strcat(handles.pathname,filesep,handles.filenametrace_unfilt));
        if fid == -1
                disp(strcat(handles.filenametrace_unfilt,' could not be loaded!!'))
        else
            handles.rawunfilteredtrace_vect = cell2mat(textscan(fid, '%f')); % read in rawtrace file
            if get(handles.check_flip,'Value')
                handles.rawunfilteredtrace_vect = -1.*handles.rawunfilteredtrace_vect;
            end
            fclose(fid); % close the file
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    filenameextrainfo = strcat(handles.pathname, filesep,'analysis_extra_info.txt');
    if isequal(exist(filenameextrainfo),2) % if the extra info file exists (legacy)
        fid = fopen(filenameextrainfo); % open it
        if fid == -1
            disp(strcat(filenameextrainfo,' could not be loaded!!'))
        else
            handles.event_more_info = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f')); % get the extra info
                % this contains :
                % 1 file number
                % 2 resistance
                % 3 time - local
                % 4 time - global
                % 5 number of local min
                % 6 index first local min
                % 7 index last local min
                % 8 dwellold (between start and stop)
                % 9 amplitude good
            fclose(fid); % close the file
        end

        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.edit_status,'String','Opening events file...done...opening rawtrace file...done')
    guidata(gcbo, handles); % update handles
    
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function listbox_eventsbox_Callback(~, ~, handles)

    set(handles.edit_status,'String','Event selected...plotting...')
    fun_plot_event('X', handles) % plot the event that was clicked
    set(handles.edit_status,'String','Event selected...plotting...done')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function listbox_eventsbox_CreateFcn(hObject, ~, ~)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_detect_local_Callback(~, ~, handles)

    set(handles.button_open,'BackgroundColor',[1 0 0])
    set(handles.edit_status,'String','Detecting local structures...')
    drawnow
    fun_find_local_struct(-1, 1, handles)
    % handles=guidata(gcbo); % update local handles
    set(handles.button_open,'BackgroundColor',[0.925 0.914 0.847])
    set(handles.edit_status,'String','Detecting local structures... done')
    fun_update_peak_stats
    % guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_remove_event_Callback(~, ~, handles)

    fun_delete_event(handles)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_delete_event(handles)
    handles.eventselected = get(handles.listbox_eventsbox,'Value'); % find event selected
    eventid = handles.event_local_struct_info_listbox(handles.eventselected,1);
    last_event = size(handles.event_local_struct_info_listbox,1);
    
    % remove event from handles.event_local_struct_info
%     index_to_delete = handles.event_local_struct_info(:,1) == eventid;
%     handles.event_local_struct_info(index_to_delete,:) = [];
%     
%     if size(handles.peak_info,1) > 0
%         % remove all corresponding peaks from handles.peak_info
%         index_to_delete = handles.peak_info(:,1) == eventid;
%         handles.peak_info(index_to_delete,:) = [];
%     end

    handles.event_local_struct_info(eventid,9) = -1;
    
    if handles.eventselected == last_event && handles.eventselected ~= 1
        
        set(handles.listbox_eventsbox,'Value',handles.eventselected - 1);
        
    end
        
    guidata(gcbo, handles); % save handles
    
    fun_update_listbox()
    
    fun_plot_event('X', handles) % plot the event that was clicked
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_method_Callback(~, ~, ~)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_method_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_baseline_Callback(~, ~, ~)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_baseline_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_status_Callback(~, ~, ~)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_status_CreateFcn(hObject, ~, ~)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_event(location, handles) % plot events in the GUI

    handles.eventselected = get(handles.listbox_eventsbox,'Value'); % find event selected

    eventid = handles.event_local_struct_info_listbox(handles.eventselected,1);

    eventpoints = [(handles.event_info(eventid,8)) (handles.event_info(eventid,9))]; % start stop points
    handles.extrapoints = handles.readextrapoints; % number of extrapoints
    
    if get(handles.check_show_extrapoints,'Value') % extra points
        
        eventtrace = handles.rawtrace_vect(eventpoints(1)-handles.extrapoints:eventpoints(2)+handles.extrapoints) - handles.event_info(eventid,5); % build event trace
        timetemp = transpose(-handles.extrapoints*handles.timestep*1000:(handles.timestep*1000):handles.timestep*(eventpoints(2)-eventpoints(1)+handles.extrapoints)*1000); % build event time vector
    
    else % no extra points
        
        eventtrace = handles.rawtrace_vect(eventpoints(1):eventpoints(2)) - handles.event_info(eventid,5); % build event trace
        timetemp = transpose(0:(handles.timestep*1000):handles.timestep*(eventpoints(2)-eventpoints(1))*1000); % build event time vector
        
    end
    
%     if get(handles.check_flip,'Value')
%     	eventtrace = eventtrace.*-1;
%     end   
    
    if strcmp(location, 'X') % plot it here
        cla(handles.axes_eventplot,'reset')
        axes(handles.axes_eventplot)
    elseif strcmp(location, 'n')
        figure() % plot in new window
    end
    baseline_DNA = str2num(get(handles.edit_baseline,'String')); % use baseline * threshold
    baseline = zeros(length(eventtrace)); % in the conductance case it is zero
    
%     plot(timetemp, baseline,'-y') % plot baseline
    
    pp4 = plot(timetemp,eventtrace,'-b'); % plot event trace
    hold on
    xlim([min(timetemp) max(timetemp)]); 
    
    if strcmp(location, 'X') % plot it here
        ylimits = get(handles.axes_eventplot,'YLim');
    else
        ylimits = get(gca,'YLim');
    end
    
    ylim([ylimits(1) ylimits(2)]);  % static limits
    
    if get(handles.check_show_region_analysis,'Value') & isfield(handles,'region_analysis_trace')
        
        if get(handles.check_show_extrapoints,'Value') % extra points
        
            eventtrace_region = handles.region_analysis_trace(eventpoints(1)-handles.extrapoints:eventpoints(2)+handles.extrapoints); % build event trace
            
        else % no extra points

            eventtrace_region = handles.region_analysis_trace(eventpoints(1):eventpoints(2)); % build event trace
            
        end
        
        plot(timetemp, eventtrace_region,'m') % plot region analysis
        
        if get(handles.check_show_thresholds,'Value')
        
            plot(timetemp, -1*baseline_DNA*str2num(get(handles.edit_region_threshold,'String')),'-c') % plot regionthreshold

        end
        
    end
    if get(handles.check_show_thresholds,'Value')
        
        plot(timetemp, -1*baseline_DNA*str2num(get(handles.edit_threshold,'String')),'-y') % plot threshold
        
    end
    
%     ylimits = get(handles.axes_eventplot,'YLim');
%     xlimits = get(handles.axes_eventplot,'XLim');
%     axis([min(timetemp) max(timetemp) ylimits(1) ylimits(2)]) % scale to fit in X

    hXLabel = xlabel('Time (ms)');
    hYLabel = ylabel('Current (nA)');
    
    % style
    set( gca                       , ...
        'FontName'   , 'Arial' );
    set([hXLabel, hYLabel], ...
        'FontName'   , 'Arial');
    set([hXLabel, hYLabel]  , ...
        'FontSize'   , 10          );
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'on'      , ...
      'YMinorTick'  , 'on'      , ...
      'YGrid'       , 'on'      , ...
      'XColor'      , [0 0 0], ...
      'YColor'      , [0 0 0], ...
      'LineWidth'   , 1         );
    set(pp4                     , ...
      'Color'           , [0 0 1], ...
      'LineWidth'       , 2           );
    

    peak_plot_info = fun_find_local_struct(eventid, 1, handles);
    
    
    plot(timetemp, -1*baseline_DNA,'-g') % plot baseline
    plot(timetemp, -2*baseline_DNA,'-g') % plot 2nd baseline
    plot(timetemp, -3*baseline_DNA,'-g') % plot 3rd baseline
    plot(timetemp, -4*baseline_DNA,'-g') % plot 4th baseline
    plot(timetemp, -5*baseline_DNA,'-g') % plot 5th baseline
    plot(timetemp, -6*baseline_DNA,'-g') % plot 6th baseline
    plot(timetemp, -7*baseline_DNA,'-g') % plot 7th baseline
    plot(timetemp, -8*baseline_DNA,'-g') % plot 8th baseline
    plot(timetemp, -9*baseline_DNA,'-g') % plot 9th baseline
    plot(timetemp, -10*baseline_DNA,'-g') % plot 10th baseline
    plot(timetemp, -11*baseline_DNA,'-g') % plot 11th baseline
    plot(timetemp, -12*baseline_DNA,'-g') % plot 12th baseline
    plot(timetemp, -13*baseline_DNA,'-g') % plot 13th baseline
    plot(timetemp, -14*baseline_DNA,'-g') % plot 14th baseline
    plot(timetemp, -15*baseline_DNA,'-g') % plot 15th baseline
%     axis([min(timetemp) max(timetemp) ylimits(1) ylimits(2)]) % scale to fit in X
    
    number_of_peaks = size(peak_plot_info,1);
    if number_of_peaks > 0
        
        if get(handles.check_list_peak_center_time,'Value') % put center peak time into box
            peak_center_time = peak_plot_info(:,4); % peak_plot_info(:,3)./2 + 
            set(handles.edit_peakdwelldata,'String',num2str(peak_center_time))
            
            clipboard('copy', num2str(peak_center_time))
            
        else % put peak dwell time into box
%             set(handles.edit_peakdwelldata,'String',num2str([peak_plot_info(:,3);peak_plot_info(:,4)]))
            set(handles.edit_peakdwelldata,'String',num2str(peak_plot_info(:,3)))
%             clipboard('copy', num2str(peak_plot_info(:,3)))
%             clipboard('copy', [num2str(peak_plot_info(1,3)),java.lang.System.getProperty('line.separator').char,num2str(peak_plot_info(1,4))])
        end
        
        
        % put giant red dots at the peaks
        if get(handles.check_show_extrapoints,'Value') % extra points
            plot(timetemp(peak_plot_info(:,2)+handles.extrapoints), eventtrace(peak_plot_info(:,2)+handles.extrapoints),'.r','MarkerSize',10) % plot red dots
        else % no extra points
            plot(timetemp(peak_plot_info(:,2)), eventtrace(peak_plot_info(:,2)+handles.extrapoints),'.r','MarkerSize',10) % plot red dots
        end
            
        %draw FWHM for each one
        for i=1:number_of_peaks
            
            plot([peak_plot_info(i,4) peak_plot_info(i,6)], [peak_plot_info(i,5) peak_plot_info(i,7)],'-m') % plot FWHM
            
        end
    end
    if strcmp(location, 's')
        dir_str = strcat(handles.pathname, filesep,'events', filesep); % in the plots dir
        if(~exist(dir_str, 'dir'))
            mkdir(dir_str);
        end
        file_str = strcat(dir_str,'event_',num2str(eventid),'_filtered.dat'); % create dat file with the data
        dlmwrite(file_str, eventtrace, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
%         if get(handles.check_overlay_unfiltered,'Value') && isfield(handles,'rawunfilteredtrace_vect')
%             file_str2 = strcat(dir_str,'event_',num2str(eventid),'_unfiltered.dat'); % create dat file with the data
%             dlmwrite(file_str2, eventtraceunfilt, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
%         end

    end
    
    % this was used in Fig 4 of the paper explaining the analysis method 
    % for Plesa et al. 2015
    extra_stuff_on = 0;
    if extra_stuff_on == 1
        
        start_end_index = find(eventtrace < -1*baseline_DNA);
        length_rect = timetemp(start_end_index(end)) - timetemp(start_end_index(1));
        xvert = [timetemp(start_end_index(1)) timetemp(start_end_index(1)) timetemp(start_end_index(end)+1) timetemp(start_end_index(end)+1)];
        yvert = [-1*baseline_DNA -2*baseline_DNA -2*baseline_DNA -1*baseline_DNA];
        ccolor = [1 1 1 1];
        p = patch(xvert,yvert,ccolor,'EdgeColor','none','FaceColor','g');
%         set(p,'FaceAlpha',0.5);

        index_high = find(eventtrace > -1*baseline_DNA & timetemp > timetemp(start_end_index(1)-1) & timetemp < timetemp(start_end_index(end)+1));
        eventtrace_mod = eventtrace;
        eventtrace_mod(index_high) = -1*baseline_DNA;
        index_low = find(eventtrace < -2*baseline_DNA & timetemp > timetemp(start_end_index(1)-1) & timetemp < timetemp(start_end_index(end)+1));
        eventtrace_mod(index_low) = -2*baseline_DNA;

        x_fill=timetemp(start_end_index(1):start_end_index(end));                  %#initialize x array
        y1_fill=-1*baseline_DNA*ones(size(x_fill));                      %#create first curve
        y2_fill=eventtrace_mod(start_end_index(1):start_end_index(end));                   %#create second curve
        y2_fill(1) = y1_fill(1);
        y2_fill(end) = y1_fill(end);
        [fillhandle,msg]=jbfill(x_fill,y1_fill,y2_fill,'red','none',1,1);
        % http://www.mathworks.nl/matlabcentral/fileexchange/13188-shade-area-between-two-curves
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function peak_plot_info = fun_find_local_struct(eventid, type_analysis, handles)
    
    % type_analysis: 1 - normal, 0 - region analysis
    baseline_DNA = str2num(get(handles.edit_baseline,'String')); % use baseline * threshold
    extrema = -1; % peaks going down
    
    peak_FWHM_interpolation_point = str2num(get(handles.edit_FWHM_level,'String')); % use 0.5 for clear peaks, use 0.3 for blockade with plateau
    
    if type_analysis == 1 % normal
        
        thresh = extrema.*str2num(get(handles.edit_threshold,'String')).*baseline_DNA; % use baseline * threshold
        selection_factor = str2num(get(handles.edit_selection,'String'));
    else % region analysis
        
        thresh = extrema.*str2num(get(handles.edit_region_threshold,'String')).*baseline_DNA; % use baseline * threshold
        selection_factor = str2num(get(handles.edit_region_selection,'String'));
    end
    
    if eventid == -1 % scan through all events
        start_for = 1;
        end_for = size(handles.event_local_struct_info,1);
        handles.peak_info = []; % reset peak info
    else % just do one event
        start_for = find(handles.event_local_struct_info(:,1)==eventid);
        end_for = start_for; 
    end
    
    peak_counter = 1;
    
    for i=start_for:end_for

        eventpoints = [(handles.event_info(handles.event_local_struct_info(i,1),8)) (handles.event_info(handles.event_local_struct_info(i,1),9))]; % start stop points
        
        eventtrace_non_region = handles.rawtrace_vect(eventpoints(1):eventpoints(2)) - handles.event_info(handles.event_local_struct_info(i,1),5); % build event trace
        if type_analysis == 1
            eventtrace = eventtrace_non_region; % build event trace
        else
            eventtrace = handles.region_analysis_trace(eventpoints(1):eventpoints(2)); % build event trace
        end
        
%         if get(handles.check_flip,'Value')
%             eventtrace = eventtrace.*-1;
%         end

        if size(eventtrace,2) > size(eventtrace,1)
            eventtrace = transpose(eventtrace);
        end
        
        event_time_vect = transpose(0:(handles.timestep*1000):handles.timestep*(eventpoints(2)-eventpoints(1))*1000); % build event time vector

        sel = (max(1.*eventtrace)-min(1.*eventtrace))/4.*selection_factor;

        [peakLoc] = peakfinder(eventtrace,sel,thresh,extrema);
%         disp(strcat('threshold: ',num2str(thresh)))
        num_peaks = size(peakLoc,1);
%       assignin('base', 'num_peaks', num_peaks)
%        assignin('base', 'eventtrace', eventtrace)
%       assignin('base', 'sel', sel)
%       assignin('base', 'thresh', thresh)
%       
%        assignin('base', 'peakLoc', peakLoc)
        
        if type_analysis == 1
            handles.event_local_struct_info(i,4) = num_peaks;
        else
            handles.event_local_struct_info(i,11) = num_peaks;
        end
            
        if num_peaks > 0 && ~isempty(peakLoc)
%             disp('here')
            % calc start and on FWHM of this event based on the baseline
            
            if type_analysis == 1
                event_trace_norm = (eventtrace) / baseline_DNA; % for FWHM calc
            else
                event_trace_norm = (eventtrace_non_region) / baseline_DNA; % for FWHM calc
            end
            
            % detect starting FWHM point
            N = length(event_trace_norm);
            icount = 2; % start at second point and go forward
            while sign(event_trace_norm(icount)+0.5) == sign(event_trace_norm(icount-1)+0.5) % go forward until you find the crossing
                icount = icount + 1;
                if icount == N + 1 %% if you cant find changeover
                    icount = N; % keep calm, but FreakOut!!
                    disp(strcat('problem with finding starting FWHM point of event #',{' '},num2str(handles.event_local_struct_info(i,1))))
                    break % and break the loop
                end
            end

            % interpolate the time of the starting FWHM point
            yvalues = [icount-1 icount];
            xvalues = [event_trace_norm(icount-1) event_trace_norm(icount)];
            tlead_index = interp1(xvalues,yvalues,-0.5);
            tlead = tlead_index*handles.timestep*1000 - handles.timestep*1000;
            
            icountend = N;   %start search for next crossing at the end
            % work your way backwards
            while ((icountend >= icount) && (sign(event_trace_norm(icountend)+0.5) == sign(event_trace_norm(icountend-1)+0.5)))
                icountend = icountend - 1;
                if icountend == 0 %% if you cant find changeover
                    icountend = 1; % keep calm, but FreakOut!!
                    disp(strcat('problem with finding ending FWHM point of event #',{' '},num2str(handles.event_local_struct_info(i,1))))
                    break % and break the loop
                end
            end
            if icountend > icount % sanity check
                
                yvalues = [icountend-1 icountend];
                xvalues = [event_trace_norm(icountend-1) event_trace_norm(icountend)];
                trail_index = interp1(xvalues,yvalues,-0.5);
                ttrail = trail_index*handles.timestep*1000 - handles.timestep*1000;
                
                event_FWHM = ttrail-tlead;% Full Width at Half Maximum in ms

            else % big problem caught here
                disp('Step-Like Pulse, no second edge, if you see this, be afraid!! Very afraid!!')
                disp(strcat('Event ID:',num2str(handles.event_local_struct_info(i,1))))
                event_FWHM = 0; % Bad case

            end
            
            % update the values just calculated
            handles.event_local_struct_info(i,3) = event_FWHM;
            handles.event_local_struct_info(i,5) = tlead;
            handles.event_local_struct_info(i,6) = ttrail;
            
            for j=1:num_peaks

                if type_analysis == 1
                    peakMag = abs(eventtrace(peakLoc(j))) - baseline_DNA;
                    if get(handles.check_basline_in_FWHM,'Value')
                        event_trace_norm = (eventtrace)./peakMag;
                    else
                        event_trace_norm = (eventtrace + baseline_DNA)./peakMag;
                    end
                    
                else
                    peakMag = abs(eventtrace(peakLoc(j))); % baseline DNA already removed
                    event_trace_norm = eventtrace./peakMag;
                end

                % detect starting FWHM point
                N = length(event_trace_norm);
                icount = peakLoc(j); % start at peak and go back
                while sign(event_trace_norm(icount)+peak_FWHM_interpolation_point) == sign(event_trace_norm(icount-1)+peak_FWHM_interpolation_point) % go backward until you find the crossing
                    icount = icount - 1;
                    if icount == 1 %% if you cant find changeover
                        icount = peakLoc(j)-1; % keep calm, but FreakOut!!
                        %disp('problem with finding starting FWHM point of event #',{' '},num2str(handles.event_local_struct_info(i,1)))
                        break % and break the loop
                    end
                end

                % crossing is between icount and icount-1
                % flip x-y axes
                % interpolate
                % flip back

                yvalues = [icount-1 icount];
                xvalues = [event_trace_norm(icount-1) event_trace_norm(icount)];
                
                timeinterp_start = interp1(xvalues,yvalues,-peak_FWHM_interpolation_point)*handles.timestep*1000 - handles.timestep*1000; 
                valueinterp_start = interp1(event_time_vect,eventtrace,timeinterp_start); 
                
                icountend = peakLoc(j);   %start search for next crossing at the end

                % work your way forward
                while sign(event_trace_norm(icountend)+peak_FWHM_interpolation_point) == sign(event_trace_norm(icountend+1)+peak_FWHM_interpolation_point)
                    icountend = icountend + 1;
                    if icountend == N %% if you cant find changeover
                        icountend = peakLoc(j)+1; % keep calm, but FreakOut!!
                        err_msg = strcat('problem with finding ending FWHM point of event ID#',{' '},num2str(handles.event_local_struct_info(i,1)));
                        disp(err_msg)
                        break % and break the loop
                    end
                end
                % cosssing is between icountend and icountend+1
                if icountend > icount % sanity check
                    
                    yvalues = [icountend icountend+1];
                    xvalues = [event_trace_norm(icountend) event_trace_norm(icountend+1)];

                    timeinterp_end = interp1(xvalues,yvalues,-peak_FWHM_interpolation_point)*handles.timestep*1000 - handles.timestep*1000; 
                    valueinterp_end = interp1(event_time_vect,eventtrace,timeinterp_end);

                    peakFWHM = timeinterp_end-timeinterp_start;% Full Width at Half Maximum in ms


                else % big problem caught here
                    
                    disp('Step-Like Pulse, no second edge, if you see this, be afraid!! Very afraid!!')
                    disp(strcat('Event ID:',num2str(handles.event_local_struct_info(i,1))))
                    disp(strcat('Peak #:',num2str(j)))
                    timeinterp_end = N*handles.timestep*1000 - handles.timestep*1000;
                    valueinterp_end = 0;
                    peakFWHM = 0; % Bad case

                end
                
                if eventid == -1
                    
                    peakpos_relative_time = (peakLoc(j)*handles.timestep*1000 - handles.timestep*1000)-tlead;

                    if get(handles.check_fwhm_center_peak_pos,'Value')
                        peakpos_relative_time = timeinterp_start;
                    end
                    
                    peakpos_relative_time_normalized = peakpos_relative_time/event_FWHM;
                    
                    if type_analysis == 1
                                                
                        
                        handles.peak_info(peak_counter,:) = [handles.event_local_struct_info(i,1) ...
                                                             peak_counter peakLoc(j) peakpos_relative_time peakFWHM ...
                                                             timeinterp_start valueinterp_start peakMag valueinterp_end ...
                                                             peakpos_relative_time_normalized];
                        % 1 - eventid
                        % 2 - peakid
                        % 3 - peakpos_ind
                        % 4 - peakpos_rel to FWHM start in time
                        % 5 - peak FWHM   
                        % 6 - peak FWHM start time
                        % 7 - peak FWHM start mag
                        % 8 - peak amp nA   
                        % 9 - peak FWHM end mag
                        % 10 - peakpos normalized with event FWHM
                    else
                        
                        
                        
                        handles.region_analysis_peak_info(peak_counter,:) = [handles.event_local_struct_info(i,1) peak_counter peakLoc(j) peakpos_relative_time peakFWHM timeinterp_start valueinterp_start peakMag valueinterp_end peakpos_relative_time_normalized];
                        % 1 - eventid
                        % 2 - peakid
                        % 3 - peakpos_ind
                        % 4 - peakpos_rel to FWHM start in time
                        % 5 - peak FWHM   
                        % 6 - peak FWHM start time
                        % 7 - peak FWHM start mag
                        % 8 - peak amp nA   
                        % 9 - peak FWHM end mag
                        % 10 - peakpos normalized with event FWHM
                        
                        % this was used to check idea
%                         handles.region_analysis_peak_info_base_cross_dwell(peak_counter,:) = peak_1st_DNA_dwell_bas;
                        
                    end
                    
                else
                    peak_plot_info(peak_counter,:) = [peak_counter peakLoc(j) peakFWHM timeinterp_start valueinterp_start timeinterp_end valueinterp_end];
                    % 1 - peakid
                    % 2 - peakpos_ind
                    % 3 - peak FWHM   
                    % 4 - peak FWHM start time
                    % 5 - peak FWHM start mag
                    % 6 - peak FWHM end time
                    % 7 - peak FWHM end mag
                end
                peak_counter = peak_counter + 1;
            end
            if type_analysis == 0 
                handles.event_local_struct_info(handles.event_local_struct_info(i,1),8) = max(handles.region_analysis_peak_info(peak_counter-num_peaks:peak_counter-1,5));
            end
        else
            if type_analysis == 0 
                handles.event_local_struct_info(handles.event_local_struct_info(i,1),8) = 0;
            end
        end
    end
    if eventid == -1
        peak_plot_info = 1;
    end
    if eventid ~= -1 && ~exist('peak_plot_info')
        peak_plot_info = [];
    end
    
    guidata(gcbo, handles); % save handles
    
    if type_analysis == 1
        fun_update_listbox()
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_threshold_Callback(~, ~, ~)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_threshold_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_selection_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_selection_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_flip_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_save_ana_Callback(~, ~, handles)

    set(handles.edit_status,'String','Saving analysis....')

    formatOut = 'YYYY_mm_dd_-_HH_MM';
    
    if ~strcmp(handles.pathname(end), filesep)
        handles.pathname = strcat(handles.pathname,filesep);
    end
    
    dir_str = strcat(handles.pathname,datestr(now,formatOut), filesep); % in this sub directory
    if(~exist(dir_str, 'dir'))
        mkdir(dir_str);
    end
    
    % copy ~detect.par
    filename = '~detect.par';
    oldname = strcat(handles.pathname, filesep, filename);
    newname = strcat(dir_str, filename);
    disp(oldname)
    disp(newname)
    copyfile(oldname,newname); % move the file
    
    % copy analysis_events
    filename = 'analysis_events';
    oldname = strcat(handles.pathname, filesep, filename);
    newname = strcat(dir_str, filename);
    copyfile(oldname,newname); % move the file
    
    % copy analysis_unfiltered_rawtrace
    filename = 'analysis_unfiltered_rawtrace';
    oldname = strcat(handles.pathname, filesep, filename);
    if isequal(exist(oldname),2)
        newname = strcat(dir_str, filename);
        copyfile(oldname,newname); % move the file
    end
    
    
    % copy analysis_rawtrace
    filename = 'analysis_rawtrace';
    oldname = strcat(handles.pathname, filesep, filename);
    newname = strcat(dir_str, filename);
    copyfile(oldname,newname); % move the file
    
    % copy analysis_parameters.txt
    filename = 'analysis_parameters.txt';
    oldname = strcat(handles.pathname, filesep, filename);
    newname = strcat(dir_str, filename);
    copyfile(oldname,newname); % move the file
    
    % copy analysis_extra_info.txt
    filename = 'analysis_extra_info.txt';
    oldname = strcat(handles.pathname, filesep, filename);
    newname = strcat(dir_str, filename);
    copyfile(oldname,newname); % move the file
    
    % save event_local_struct_info
    file_str = strcat(dir_str,'event_local_struct_info');
    dlmwrite(file_str, handles.event_local_struct_info, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    
    % save peak_info
    file_str = strcat(dir_str,'peak_info');
    dlmwrite(file_str, handles.peak_info, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    
    % save region_analysis_peak_info
    file_str = strcat(dir_str,'region_analysis_peak_info');
    dlmwrite(file_str, handles.region_analysis_peak_info, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
       
    fun_save_parameters(dir_str, handles)
    set(handles.edit_status,'String','Saving analysis.... done')
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_save_parameters(location, handles)

    if location == -1 % ask where to save
        savepath = uigetdir(handles.pathname, 'Pick a Directory');
    else %
        savepath = location;
    end

    file_str = strcat(savepath, filesep,'~local_structure_analysis.par');
    fid = fopen(file_str, 'w');
    fprintf(fid, '%s\r\n', '%% Local Structure Analysis Par File');
    fprintf(fid, '%s\r\n', strcat('%% Version 0.4 - March 2013'));
    fprintf(fid, '%s\r\n', '%% Flip current');
    fprintf(fid, '%s\r\n', num2str(get(handles.check_flip,'Value')));
    fprintf(fid, '%s\r\n', '%% Method');
    fprintf(fid, '%s\r\n', num2str(get(handles.pop_method,'Value')));
    fprintf(fid, '%s\r\n', '%% Baseline - nA');
    fprintf(fid, '%s\r\n', get(handles.edit_baseline,'String'));
    fprintf(fid, '%s\r\n', '%% PeakFinder Threshold');
    fprintf(fid, '%s\r\n', get(handles.edit_threshold,'String'));
    fprintf(fid, '%s\r\n', '%% PeakFinder Selection');
    fprintf(fid, '%s\r\n', get(handles.edit_selection,'String'));
    
    fprintf(fid, '%s\r\n', '%% Region Analysis - paramA');
    fprintf(fid, '%s\r\n', get(handles.edit_region_param_a,'String'));
    fprintf(fid, '%s\r\n', '%% Region Analysis - Threshold');
    fprintf(fid, '%s\r\n', get(handles.edit_region_threshold,'String'));
    fprintf(fid, '%s\r\n', '%% Region Analysis - Selection');
    fprintf(fid, '%s\r\n', get(handles.edit_region_selection,'String'));
    
    fprintf(fid, '%s\r\n', '%% Dwell Bins');
    fprintf(fid, '%s\r\n', get(handles.edit_dwell_bins,'String'));
    fprintf(fid, '%s\r\n', '%% Amp Bins');
    fprintf(fid, '%s\r\n', get(handles.edit_amp_bins,'String'));
    fprintf(fid, '%s\r\n', '%% Pos Bins');
    fprintf(fid, '%s\r\n', get(handles.edit_pos_bins,'String'));
    fprintf(fid, '%s\r\n', '%% Current Bins');
    fprintf(fid, '%s\r\n', get(handles.edit_current_bins,'String'));
    fprintf(fid, '%s\r\n', '%% EOF');
    fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_load_parameters(location, handles)

    if location == -1 % ask where from
        [FileName,PathName] = uigetfile('*.par','Select the MATLAB code file'); % get PAR location from user
        loadpath = strcat(PathName, FileName); % full path and filename
    else %
        loadpath = location;
    end

    if isequal(exist(loadpath,'file'),2) % 
        fid = fopen(loadpath);
        if fid == -1
            disp(strcat(loadpath,' could not be loaded!!'))
        else
            line = fgetl(fid); %% Local Structure Analysis Par File
            line = fgetl(fid); %% Version #
            line = fgetl(fid); %% Flip current
            flipcurrent = str2double(fgetl(fid)); %%
            if flipcurrent
                set(handles.check_flip,'Value',1);
            else
                set(handles.check_flip,'Value',0);
            end
            line = fgetl(fid); %% Method
            methodtype = str2double(fgetl(fid)); %% 
            set(handles.pop_method,'Value',methodtype);
            line = fgetl(fid); %% Baseline
            set(handles.edit_baseline,'String',fgetl(fid));
            line = fgetl(fid); %% PeakFinder Threshold
            set(handles.edit_threshold,'String',fgetl(fid));
            line = fgetl(fid); %% PeakFinder Selection
            set(handles.edit_selection,'String',fgetl(fid));
            
            line = fgetl(fid); %% Region Analysis - paramA
            set(handles.edit_region_param_a,'String',fgetl(fid));
            line = fgetl(fid); %% Region Analysis - Threshold
            set(handles.edit_region_threshold,'String',fgetl(fid));
            line = fgetl(fid); %% Region Analysis - Selection
            set(handles.edit_region_selection,'String',fgetl(fid));
            
            line = fgetl(fid); %% Dwell Bins
            set(handles.edit_dwell_bins,'String',fgetl(fid));
            line = fgetl(fid); %% Amp Bins
            set(handles.edit_amp_bins,'String',fgetl(fid));
            line = fgetl(fid); %% Pos Bins
            set(handles.edit_pos_bins,'String',fgetl(fid));
            line = fgetl(fid); %% Current Bins
            set(handles.edit_current_bins,'String',fgetl(fid));
            
        end
        fclose(fid);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_load_ana_Callback(hObject, ~, ~)
    
    loadpath = uigetdir(pwd, 'Pick a Directory');
    
    fun_load_analysis_files(loadpath)
    
    handles=guidata(gcbo); % update local handles

    % read in par file
    
    fun_load_parameters(strcat(loadpath,filesep,'~local_structure_analysis.par'), handles)
    
    handles=guidata(gcbo); % update local handles

    % read in event_local_struct_info
    filenamelocalstructinfo = strcat(loadpath, filesep,'event_local_struct_info');
    if isequal(exist(filenamelocalstructinfo,'file'),2) % if the file exists 
        fid = fopen(filenamelocalstructinfo); % open it
        if fid == -1
            disp(strcat(filenamelocalstructinfo,' could not be loaded!!'))
        else
            handles.event_local_struct_info = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f')); % get the info
            % 1 id
            % 2 integral
            % 3 FWHM_dwell
            % 4 number of local struct
            % 5 event start ind
            % 6 event end ind
            % 7 region integral
            % 8 max peak FWHM
            % 9 group 1 type:
            %       -1 - bad
            %       0 - unsorted
            %       1 - folded
            %       2 - clogged
            %       3 - unfolded
            %       
            % 10 group 2 type
            %       0 - unsorted
            %       1 - bare (no spikes)
            %       2 - spike A
            %       3 - spike B
            %       4 - spike C
            %       5 - spike D
            %       6 - spike E
            %       7 - spike F
            %       8 - spike G
            %       9 - spike H
            %       10 - spike I
            % 11 number of local struct from region analysis
            % 12 mean first extra points
            % 13 mean last extra points
            % 14 selected 1 -yes 0 - no
        end
        fclose(fid);
    end
        
    % read in peak_info
    filenamepeakinfo = strcat(loadpath, filesep,'peak_info');
    if isequal(exist(filenamepeakinfo,'file'),2) % if the file exists 
        fid = fopen(filenamepeakinfo); % open it
        if fid == -1
            disp(strcat(filenamepeakinfo,' could not be loaded!!'))
        else
            handles.peak_info = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f %f')); % get the extra info
            % 1 - eventid
            % 2 - peakid
            % 3 - peakpos_ind
            % 4 - peakpos_rel to FWHM start in time
            % 5 - peak FWHM   
            % 6 - peak FWHM start time
            % 7 - peak FWHM start mag
            % 8 - peak amp nA   
            % 9 - peak FWHM end mag
            % 10 - peakpos normalized with event FWHM
        end
        fclose(fid);
    end
    
        
    % read in region_analysis_peak_info
    filenamepeakinfo = strcat(loadpath, filesep,'region_analysis_peak_info');
    if isequal(exist(filenamepeakinfo,'file'),2) % if the file exists 
        fid = fopen(filenamepeakinfo); % open it
        if fid == -1
            disp(strcat(filenamepeakinfo,' could not be loaded!!'))
        else
            handles.region_analysis_peak_info = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f %f')); % get the extra info
            % 1 - eventid
            % 2 - peakid
            % 3 - peakpos_ind
            % 4 - peakpos_rel to FWHM start in time
            % 5 - peak FWHM   
            % 6 - peak FWHM start time
            % 7 - peak FWHM start mag
            % 8 - peak amp nA   
            % 9 - peak FWHM end mag
            % 10 - peakpos normalized with event FWHM
        end
        fclose(fid);
    end
    
    guidata(hObject, handles); % update handles
    fun_update_listbox()
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_plot_individual_Callback(~, ~, handles)

    %new window code here
    set(handles.edit_status,'String','Generating New Figure');

    % find what type of plot
    selection = get(handles.uipanel2,'SelectedObject');

    % call the appropriate plotting function with the 'n' tag to indicate new
    % figure

    switch get(selection,'Tag')
        case 'radio_dwell'
            fun_plot_dwell(handles, 'n')
        case 'radio_amp'
            fun_plot_amplitude(handles, 'n')
        case 'radio_pos'
            fun_plot_pos(handles, 'n')
        case 'radio_pos_abs'
            fun_plot_pos_abs(handles, 'n')
        case 'radio_current'
            fun_plot_current(handles, 'n')
        case 'radio_peaks_per_event'
            fun_plot_ppe(handles, 'n')
        case 'radio_scatter_spike'
            fun_plot_scatter_spike(handles, 'n')
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in button_mass_export.
function button_mass_export_Callback(hObject, ~, handles)

    set(handles.edit_status,'String','Mass exporting all plots to events file directory....')

    % this plots all of the different types of plots available and saves them
    % in .fig and .png format in the \plots directory

    dir_str = strcat(handles.pathname,'plots', filesep); % in this sub directory
    if(~exist(dir_str, 'dir'))
        mkdir(dir_str);
    end

    % this is to see if export_fig is available
    expfig_on = 0;
    if exist('export_fig') == 2
        expfig_on = 1;
    end
    
    % dwell hist plot
    newFig = figure;
    fun_plot_dwell(handles, 's') % data file saved inside this function
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_dwell_hist.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_dwell_hist_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_dwell_hist_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_localstruct_dwell_hist.eps')) 
    
    if expfig_on
        set(gcf, 'Color', 'w');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_dwell_hist_expfig.eps'), '-cmyk');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_dwell_hist_expfig.png'), '-cmyk');
    end
    
    close(newFig)

    
    
    % amplitude histogram plot
    newFig = figure;
    fun_plot_amplitude(handles, 's') % data file saved inside this function
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_amp_hist.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_amp_hist_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_amp_hist_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_localstruct_amp.eps')) 
    if expfig_on
        set(gcf, 'Color', 'w');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_amp_expfig.eps'), '-cmyk');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_amp_expfig.png'), '-cmyk');
    end
    close(newFig)

%     % current hist plot
%     if isfield(handles,'rawtrace_vect') % if rawtrace is loaded
%         newFig = figure;
%         fun_plot_current(handles, 's') % data file saved inside this function
%         saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_current_hist.fig'), 'fig'); % save as .fig
%         saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_current_hist_large.png'), 'png');
%         set(newFig, 'PaperPositionMode', 'auto');
%         saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_current_hist_small.png'), 'png');
%         print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_localstruct_current_hist.eps')) 
%         if expfig_on
%             set(gcf, 'Color', 'w');
%             export_fig(newFig, 'plot_localstruct_current_hist_expfig.eps', '-cmyk');
%             export_fig(newFig, 'plot_localstruct_current_hist_expfig.png', '-cmyk');
%         end
%         close(newFig)
%     end

    % plot pos hist
    newFig = figure;
    fun_plot_pos(handles, 's')
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_pos_hist.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_pos_hist_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_pos_hist_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_localstruct_pos_hist.eps'))
    if expfig_on
        set(gcf, 'Color', 'w');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_pos_hist_expfig.eps'), '-cmyk');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_pos_hist_expfig.png'), '-cmyk');
    end
    close(newFig)
    
    % plot pos hist abs
    newFig = figure;
    fun_plot_pos_abs(handles, 's')
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_pos_abs_hist.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_pos_abs_hist_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_pos_abs_hist_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_localstruct_pos_abs_hist.eps'))
    if expfig_on
        set(gcf, 'Color', 'w');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_pos_abs_hist_expfig.eps'), '-cmyk');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_pos_abs_hist_expfig.png'), '-cmyk');
    end
    close(newFig)
    
    % plot ppe hist
    newFig = figure;
    fun_plot_ppe(handles, 's')
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_ppe_hist.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_ppe_hist_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_ppe_hist_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_localstruct_ppe_hist.eps'))
    if expfig_on
        set(gcf, 'Color', 'w');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_ppe_hist_expfig.eps'), '-cmyk');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_ppe_hist_expfig.png'), '-cmyk');
    end
    close(newFig)

    % plot pop dist hist
    newFig = figure;
    fun_plot_pop_dist(handles, 's')
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_pop_dist_pie.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_pop_dist_pie_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_pop_dist_pie_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_localstruct_pop_dist_pie.eps'))
    if expfig_on
        set(gcf, 'Color', 'w');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_pop_dist_pie_expfig.eps'), '-cmyk');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_pop_dist_pie_expfig.png'), '-cmyk');
    end
    close(newFig)
    
    % plot time hist
    newFig = figure;
    fun_plot_time(handles, 's')
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_time.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_time_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_time_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_localstruct_time.eps'))
    if expfig_on
        set(gcf, 'Color', 'w');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_time_expfig.eps'), '-cmyk');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_time_expfig.png'), '-cmyk');
    end
    close(newFig)
    
    % plot scatter spike
    newFig = figure;
    fun_plot_scatter_spike(handles, 's')
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_scatter_spike.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_scatter_spike_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_scatter_spike_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_localstruct_scatter_spike.eps'))
    if expfig_on
        set(gcf, 'Color', 'w');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_scatter_spike_expfig.eps'), '-cmyk');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_localstruct_scatter_spike_expfig.png'), '-cmyk');
    end
    close(newFig)
    
    % all done here

    set(handles.edit_status,'String','Mass exporting all plots to events file directory....done')

    guidata(hObject, handles); % update handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_dwell_bins_Callback(~, ~, ~)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_dwell_bins_CreateFcn(hObject, ~, ~)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_amp_bins_Callback(~, ~, ~)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_amp_bins_CreateFcn(hObject, ~, ~)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pos_bins_Callback(~, ~, ~)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pos_bins_CreateFcn(hObject, ~, ~)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_total_events_Callback(~, ~, ~)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_total_events_CreateFcn(hObject, ~, ~)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_event_plot_Callback(~, ~, handles)

    fun_plot_event('n', handles) % open it in new figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_current_bins_Callback(~, ~, ~)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_current_bins_CreateFcn(hObject, ~, ~)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_current(handles, location) % current histogram function

    fun_gen_selected_spikes_list()
        handles = guidata(gcbo);

    if handles.eventid_index_selected ~= -1
        if isfield(handles,'rawtrace_vect') % only run if rawtrace file is loaded
            colormap('default')
            if strcmp(location,'X') % plot in GUI
                cla(handles.axes_graphs,'reset')
                axes(handles.axes_graphs)
            end
            if strcmp(location,'n') % make new figure
                ylimits = get(handles.axes_graphs,'YLim');
                xlimits = get(handles.axes_graphs,'XLim');
                figure()
            end

            bins = str2num(get(handles.edit_current_bins,'String')); % number of bins

            %%%%%%%%%%%%%%%%%%% TODO, use only events with peaks
            fun_make_new_trace(handles) % reformat trace so that baselines are zero and cutoff extrapoints if necessary
            handles=guidata(gcbo); % update local handles
            %%%%%%%%%%%%%%%%%%%

        %     if get(handles.check_conductance,'Value') % convert to nS
        %     	voltage = str2double(get(handles.edit_voltage,'String'))/1000;
        %     	handles.newtrace = handles.newtrace./voltage;
        %     end


            [n,xout] = hist(handles.newtrace,bins); % make the histogram
            hist(handles.newtrace,bins)
            binwidth = xout(2)-xout(1);
        %     set(handles.edit_2_bin_width,'String',num2str(binwidth))

            % next four lines are needed to make the histogram show properly under
            % log scale
            ph = get(gca,'children');
            vn = get(ph,'Vertices');
            vn(:,2) = vn(:,2) + 1;
            set(ph,'Vertices',vn);

            set(gca,'YScale','log') % use log scale on Y

        %     if get(handles.check_conductance,'Value') % use right X label
        %         hXLabel = xlabel('Conductance (nS)');
        %     else
                hXLabel = xlabel('Current (nA)');
        %     end

            hYLabel = ylabel('Counts'); % Y label

            % style
            set( gca                       , ...
                'FontName'   , 'Arial' );
            set([hXLabel, hYLabel], ...
                'FontName'   , 'Arial');
            set([hXLabel, hYLabel]  , ...
                'FontSize'   , 10          );

            set(gca, ...
              'Box'         , 'on'     , ...
              'TickDir'     , 'in'     , ...
              'TickLength'  , [.02 .02] , ...
              'XMinorTick'  , 'on'      , ...
              'YMinorTick'  , 'on'      , ...
              'YGrid'       , 'off'      , ...
              'XColor'      , [0 0 0], ...
              'YColor'      , [0 0 0], ...
              'LineWidth'   , 1         );

            if strcmp(location,'X')
                % get limits and update limit text boxes
                xlimits = get(handles.axes_graphs,'XLim');
                ylimits = get(handles.axes_graphs,'YLim');
                set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
                set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
                set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
                set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
            end

            if strcmp(location,'n') % use old limits in new figure, for opening plot in window
                axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
            end

            if strcmp(location,'s') % mass-export case
                dir_str = strcat(handles.pathname,'plots', filesep); % in the plots directory
                if(~exist(dir_str, 'dir'))
                    mkdir(dir_str);
                end
                file_str = strcat(dir_str,'plot_data_current_raw_trace_baseline_at_zero.dat'); % make a file
                % with the new trace data
                dlmwrite(file_str, handles.newtrace, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');

                file_str = strcat(dir_str,'plot_data_current_histogram_values.dat'); % make a file
                % with the new trace histogram data
                dlmwrite(file_str, [xout' n'], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
            end
        else
            set(handles.text_status,'String','Rawtrace Not Loaded') % catch no rawtrace file case
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_make_new_trace(handles) 
    % take only the trace for selected events

    indpos = 1; % pointer to current position, initialize
    handles.extrapoints = 50; % get number of extra points to each side of the event
    handles.newtrace = 0;
    newstart = 0;

    number_of_selected_events = size(handles.event_local_struct_info,1);

    % pre allocate for speed
    newstart = zeros(number_of_selected_events);
    newend = zeros(number_of_selected_events);

    for i=1:number_of_selected_events % for number of selected events

        eventpoints = [(handles.event_info(handles.event_local_struct_info(i,1),8)) (handles.event_info(handles.event_local_struct_info(i,1),9))]; % get event start stop

        eventlength = handles.event_info(handles.event_local_struct_info(i,1),9) - handles.event_info(handles.event_local_struct_info(i,1),8) + handles.extrapoints*2 + 1; % calc length in points
        eventtrace = handles.rawtrace_vect(eventpoints(1)-handles.extrapoints:eventpoints(2)+handles.extrapoints) - handles.event_info(handles.event_local_struct_info(i,1),5); % get the event trace, and remove baseline
        newstart(i) = indpos + handles.extrapoints; % calculate the new start position

        handles.newtrace(indpos:indpos+eventlength-1) = eventtrace; % add this trace to new trace

        newend(i) = newstart(i) + handles.event_info(handles.event_local_struct_info(i,1),9) - handles.event_info(handles.event_local_struct_info(i,1),8); % calculate the new end position
        indpos = indpos + eventlength; % calc new start index value
    end
    guidata(gcbo, handles); % save handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_current_Callback(hObject, ~, handles)

    set(handles.radio_scatter_spike,'Value',0);
    set(handles.radio_dwell,'Value',0);
    set(handles.radio_amp,'Value',0);
    set(handles.radio_pos,'Value',0);
    set(handles.radio_pos_abs,'Value',0);
    set(handles.radio_current,'Value',1);
    set(handles.radio_peaks_per_event,'Value',0);
    fun_plot_current(handles,'X') % plot current hist
    guidata(hObject, handles); % update handles
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_dwell_Callback(hObject, ~, handles)

    set(handles.radio_scatter_spike,'Value',0);
    set(handles.radio_dwell,'Value',1);
    set(handles.radio_amp,'Value',0);
    set(handles.radio_pos,'Value',0);
    set(handles.radio_pos_abs,'Value',0);
    set(handles.radio_current,'Value',0);
    set(handles.radio_peaks_per_event,'Value',0);
    fun_plot_dwell(handles,'X') % plot current hist
    guidata(hObject, handles); % update handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_dwell(handles, location)


    fun_gen_selected_spikes_list()
    handles = guidata(gcbo);
    
    if handles.index_peaks_to_analyze ~= -1

        colormap('default')
        if strcmp(location,'X') % plot here
            cla(handles.axes_graphs,'reset')
            axes(handles.axes_graphs)
        end
        if strcmp(location,'n') % plot in new figure
            ylimits = get(handles.axes_graphs,'YLim');
            xlimits = get(handles.axes_graphs,'XLim');
            figure()
        end

        bins = str2num(get(handles.edit_dwell_bins,'String'));

        dwelldata = handles.peak_info(handles.index_peaks_to_analyze,5);
        pos_norm = handles.peak_info(handles.index_peaks_to_analyze,10);
        
        plot_special = 0;
        
        if plot_special
        
            dwelldata = dwelldata(pos_norm > 0.1);
            
        else
            
            
        end

        [n,xout] = hist(dwelldata,bins); % plot it
        hist(dwelldata,bins)
        binwidth = xout(2)-xout(1);
        set(handles.edit_2_bin_width,'String',num2str(binwidth))

        hXLabel = xlabel('Time (ms)');
        hYLabel = ylabel('Counts');
        
        if plot_special
        
            xlim([0 0.1])
            ylim([0 60])
            
        else
            
            
        end

        % style
        set( gca                       , ...
            'FontName'   , 'Arial' );
        set([hXLabel, hYLabel], ...
            'FontName'   , 'Arial');
        set([hXLabel, hYLabel]  , ...
            'FontSize'   , 10          );
        set(gca, ...
          'Box'         , 'on'     , ...
          'TickDir'     , 'in'     , ...
          'TickLength'  , [.02 .02] , ...
          'XMinorTick'  , 'on'      , ...
          'YMinorTick'  , 'on'      , ...
          'YGrid'       , 'on'      , ...
          'XColor'      , [0 0 0], ...
          'YColor'      , [0 0 0], ...
          'LineWidth'   , 1         );

        if strcmp(location,'X')
            % get limits and update limit text boxes
            xlimits = get(handles.axes_graphs,'XLim');
            ylimits = get(handles.axes_graphs,'YLim');
            set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
            set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
            set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
            set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
        end

        if strcmp(location,'n')% use same limits
            axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
        end

        if strcmp(location,'s')
            dir_str = strcat(handles.pathname,'plots', filesep); % in the plots directory
            if(~exist(dir_str, 'dir'))
                mkdir(dir_str);
            end
            file_str = strcat(dir_str,'plot_data_dwell.dat'); % save dat file with dwell data
            dlmwrite(file_str, dwelldata, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');

            file_str = strcat(dir_str,'plot_data_dwell_hist.dat'); % create dat file with the hist data
            dlmwrite(file_str, [xout' n'], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_pos(handles, location)


    fun_gen_selected_spikes_list()
    handles = guidata(gcbo);
    
    if handles.index_peaks_to_analyze ~= -1

        colormap('default')
        if strcmp(location,'X') % plot here
            cla(handles.axes_graphs,'reset')
            axes(handles.axes_graphs)
        end
        if strcmp(location,'n') % plot in new figure
            ylimits = get(handles.axes_graphs,'YLim');
            xlimits = get(handles.axes_graphs,'XLim');
            figure()
        end

        bins = str2num(get(handles.edit_pos_bins,'String'));

        posdata = handles.peak_info(handles.index_peaks_to_analyze,10);
        dwelldata = handles.peak_info(handles.index_peaks_to_analyze,5);
        
%         posdata = posdata(dwelldata < 0.025);

        [n,xout] = hist(posdata,bins); % plot it
        hist(posdata,bins)
        binwidth = xout(2)-xout(1);
        set(handles.edit_2_bin_width,'String',num2str(binwidth))

        hXLabel = xlabel('Normalized Position');
        hYLabel = ylabel('Counts');

        % style
        set( gca                       , ...
            'FontName'   , 'Arial' );
        set([hXLabel, hYLabel], ...
            'FontName'   , 'Arial');
        set([hXLabel, hYLabel]  , ...
            'FontSize'   , 10          );
        set(gca, ...
          'Box'         , 'on'     , ...
          'TickDir'     , 'in'     , ...
          'TickLength'  , [.02 .02] , ...
          'XMinorTick'  , 'on'      , ...
          'YMinorTick'  , 'on'      , ...
          'YGrid'       , 'on'      , ...
          'XColor'      , [0 0 0], ...
          'YColor'      , [0 0 0], ...
          'LineWidth'   , 1         );

        if strcmp(location,'X')
            % get limits and update limit text boxes
            xlimits = get(handles.axes_graphs,'XLim');
            ylimits = get(handles.axes_graphs,'YLim');
            set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
            set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
            set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
            set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
        end

        if strcmp(location,'n')% use same limits
            axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
        end

        if strcmp(location,'s')
            dir_str = strcat(handles.pathname,'plots', filesep); % in the plots directory
            if(~exist(dir_str, 'dir'))
                mkdir(dir_str);
            end
            file_str = strcat(dir_str,'plot_data_pos.dat'); % save dat file with dwell data
            dlmwrite(file_str, posdata, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');

            file_str = strcat(dir_str,'plot_data_pos_hist.dat'); % create dat file with the hist data
            dlmwrite(file_str, [xout' n'], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_2_bin_width_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_2_bin_width_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uipanel2_SelectionChangeFcn(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_amp_CreateFcn(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_current_CreateFcn(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_total_peaks_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_total_peaks_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_amp_Callback(hObject, ~, handles)

    set(handles.radio_scatter_spike,'Value',0);
    set(handles.radio_dwell,'Value',0);
    set(handles.radio_amp,'Value',1);
    set(handles.radio_pos,'Value',0);
    set(handles.radio_pos_abs,'Value',0);
    set(handles.radio_current,'Value',0);
    set(handles.radio_peaks_per_event,'Value',0);
    fun_plot_amplitude(handles,'X') % plot current hist
    guidata(hObject, handles); % update handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_amplitude(handles, location) % plot amplitude histogram


    fun_gen_selected_spikes_list()
    handles = guidata(gcbo);
    
    if handles.index_peaks_to_analyze ~= -1
    
        colormap('default')
        if strcmp(location,'X') % plot in GUI
            cla(handles.axes_graphs,'reset')
            axes(handles.axes_graphs)
        end

        if strcmp(location,'n') % plot in new figure
            ylimits = get(handles.axes_graphs,'YLim');
            xlimits = get(handles.axes_graphs,'XLim');
            figure()
        end

        bins = str2num(get(handles.edit_amp_bins,'String'));

        amplitude = handles.peak_info(handles.index_peaks_to_analyze,8); % in nA
        
        if get(handles.check_conductance,'Value')
            amplitude = amplitude./(handles.voltage/1000);
        end
        
        pos_norm = handles.peak_info(handles.index_peaks_to_analyze,10);
        
        plot_special = 0;
        
        if plot_special
            
             amplitude = amplitude(pos_norm > 0.1);
        
             xout_input = [0.182022743899998,0.203745131699998,0.225467519499998,0.247189907299998,0.268912295099998,0.290634682899998,0.312357070699998,0.334079458499998,0.355801846299998,0.377524234099998,0.399246621899998,0.420969009699998,0.442691397499998,0.464413785299998,0.486136173099998,0.507858560899998,0.529580948699998,0.551303336499998,0.573025724299998,0.594748112099998,0.616470499899998,0.638192887699998,0.659915275499998,0.681637663299998,0.703360051099998,0.725082438899998,0.746804826699998,0.768527214499998,0.790249602299998,0.811971990099998,0.833694377899998,0.855416765699998,0.877139153499998,0.898861541299998,0.920583929099998,0.942306316899998,0.964028704699998,0.985751092499998,1.00747348030000,1.02919586810000,1.05091825590000,1.07264064370000,1.09436303150000,1.11608541930000,1.13780780710000,1.15953019490000,1.18125258270000,1.20297497050000,1.22469735830000,1.24641974610000;];
             [n,xout] = hist(amplitude,xout_input); % plot it
             hist(amplitude,xout_input)  
             xlim([0 1.4])
             ylim([0 140])
             
        else
                
            [n,xout] = hist(amplitude,bins); % plot it
        %         assignin('base', 'xout_amp', xout)
        
            hist(amplitude,bins)
%             xlim([0 1.4])
%              ylim([0 140])
        end
        
        binwidth = xout(2)-xout(1);
        set(handles.edit_2_bin_width,'String',num2str(binwidth))

         
        
        if get(handles.check_conductance,'Value')
            hXLabel = xlabel('Conductance (nS)');
        else
            hXLabel = xlabel('Current Blockade (nA)');
        end

        hYLabel = ylabel('Counts');

        % style
        set( gca                       , ...
            'FontName'   , 'Arial' );
        set([hXLabel, hYLabel], ...
            'FontName'   , 'Arial');
        set([hXLabel, hYLabel]  , ...
            'FontSize'   , 10          );

        set(gca, ...
          'Box'         , 'on'     , ...
          'TickDir'     , 'in'     , ...
          'TickLength'  , [.02 .02] , ...
          'XMinorTick'  , 'on'      , ...
          'YMinorTick'  , 'on'      , ...
          'YGrid'       , 'on'      , ...
          'XColor'      , [0 0 0], ...
          'YColor'      , [0 0 0], ...
          'LineWidth'   , 1         );

        if strcmp(location,'X')
            % get limits and update limit text boxes
            xlimits = get(handles.axes_graphs,'XLim');
            ylimits = get(handles.axes_graphs,'YLim');
            set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
            set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
            set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
            set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
        end

        if strcmp(location,'n') % set new figure limits
            axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
        end

        if strcmp(location,'s')
            dir_str = strcat(handles.pathname,'plots', filesep); % in the plots dir
            if(~exist(dir_str, 'dir'))
                mkdir(dir_str);
            end
            file_str = strcat(dir_str,'plot_data_amplitude.dat'); % create dat file with the data
            dlmwrite(file_str, amplitude, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');

            file_str = strcat(dir_str,'plot_data_amplitude_hist.dat'); % create dat file with the hist data
            dlmwrite(file_str, [xout' n'], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_pos_Callback(hObject, ~, handles)

set(handles.radio_scatter_spike,'Value',0);
set(handles.radio_dwell,'Value',0);
set(handles.radio_amp,'Value',0);
set(handles.radio_pos,'Value',1);
set(handles.radio_pos_abs,'Value',0);
set(handles.radio_current,'Value',0);
set(handles.radio_peaks_per_event,'Value',0);
fun_plot_pos(handles, 'X')
handles=guidata(gcbo); % get handles
guidata(hObject, handles); % update handles
fun_update_peak_stats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_plot2xmin_Callback(hObject, ~, handles)


    % set the min X axis value

    set(handles.edit_status,'String','Setting Plot2 Xmin')
    xlimits = get(handles.axes_graphs,'XLim');
    minvalue = str2double(get(handles.edit_plot2xmin,'String'));
    set(handles.axes_graphs,'XLim',[minvalue xlimits(2)]);
    guidata(hObject, handles); % update handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_plot2xmin_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_plot2xmax_Callback(hObject, ~, handles)


    % set the max X axis value

    set(handles.edit_status,'String','Setting Plot2 Xmax')
    xlimits = get(handles.axes_graphs,'XLim');
    maxvalue = str2double(get(handles.edit_plot2xmax,'String'));
    set(handles.axes_graphs,'XLim',[xlimits(1) maxvalue]);
    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_plot2xmax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_plot2ymin_Callback(hObject, ~, handles)


    % set plot min Y value
    
    set(handles.edit_status,'String','Changing Plot2 YMin')
    ylimits = get(handles.axes_graphs,'YLim');
    minvalue = str2double(get(handles.edit_plot2ymin,'String'));
    set(handles.axes_graphs,'YLim',[minvalue ylimits(2)]);
    guidata(hObject, handles); % update handles

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_plot2ymin_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_plot2ymax_Callback(hObject, ~, handles)

    % set plot max Y value
    
    set(handles.edit_status,'String','Changing Plot2 YMax')
    ylimits = get(handles.axes_graphs,'YLim');
    maxvalue = str2double(get(handles.edit_plot2ymax,'String'));
    set(handles.axes_graphs,'YLim',[ylimits(1) maxvalue]);
    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_plot2ymax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figure1_WindowKeyReleaseFcn(hObject, eventdata, handles)
    % empty




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function listbox_eventsbox_KeyPressFcn(~, eventdata, handles)
% hObject    handle to listbox_eventsbox (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

handles.eventselected = get(handles.listbox_eventsbox,'Value'); % find event selected

eventid = handles.event_local_struct_info_listbox(handles.eventselected,1);

% assignin('base', 'key', eventdata.Key)
switch eventdata.Key
    
    case 'delete'
        fun_delete_event(handles)
        handles = guidata(gcbo);
    case 'x' % Unfolded - Spike A
        handles.event_local_struct_info(eventid,9) = 3;
        handles.event_local_struct_info(eventid,10) = 2;
    case 'm' % folded
        handles.event_local_struct_info(eventid,9) = 1;
        
    case 's' % unfolded - bare
        handles.event_local_struct_info(eventid,9) = 3;
        handles.event_local_struct_info(eventid,10) = 1;
        
end
% assignin('base', 'key', eventdata.Key)
% assignin('base', 'eventid', eventid)
guidata(gcbo, handles); % save handles
fun_update_listbox()
% handles=guidata(gcbo); % get handles
% fun_plot_event('X', handles) % plot the next event



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_par_save_Callback(~, ~, handles)

fun_save_parameters(-1, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_load_par_Callback(~, ~, handles)

fun_load_parameters(-1,handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_sort_events_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_sort_events_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_refresh_listbox_Callback(hObject, eventdata, handles)

fun_update_listbox()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_ppe_bins_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_ppe_bins_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_peaks_per_event_Callback(~, ~, handles)


set(handles.radio_scatter_spike,'Value',0);
set(handles.radio_dwell,'Value',0);
set(handles.radio_amp,'Value',0);
set(handles.radio_pos,'Value',0);
set(handles.radio_pos_abs,'Value',0);
set(handles.radio_current,'Value',0);
set(handles.radio_peaks_per_event,'Value',1);
fun_plot_ppe(handles, 'X')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_ppe(handles, location) % plot amplitude histogram


    fun_gen_selected_spikes_list()
    handles = guidata(gcbo);
    
    if handles.index_peaks_to_analyze ~= -1

        colormap('default')
        if strcmp(location,'X') % plot in GUI
            cla(handles.axes_graphs,'reset')
            axes(handles.axes_graphs)
        end

        if strcmp(location,'n') % plot in new figure
            ylimits = get(handles.axes_graphs,'YLim');
            xlimits = get(handles.axes_graphs,'XLim');
            figure()
        end

        bins = str2num(get(handles.edit_ppe_bins,'String'));
        bins_vector = 0:1:bins;


        ppe = handles.event_local_struct_info(handles.eventid_index_selected,4);


        [n,xout] = hist(ppe,bins_vector); % plot it
        hist(ppe,bins_vector)
        binwidth = xout(2)-xout(1);
        set(handles.edit_2_bin_width,'String',num2str(binwidth))

    %     if get(handles.check_conductance,'Value')
    %         hXLabel = xlabel('Conductance (nS)');
    %     else
            hXLabel = xlabel('Peaks per Event');
    %     end

        hYLabel = ylabel('Counts');

        % style
        set( gca                       , ...
            'FontName'   , 'Arial' );
        set([hXLabel, hYLabel], ...
            'FontName'   , 'Arial');
        set([hXLabel, hYLabel]  , ...
            'FontSize'   , 10          );

        set(gca, ...
          'Box'         , 'on'     , ...
          'TickDir'     , 'in'     , ...
          'TickLength'  , [.02 .02] , ...
          'XMinorTick'  , 'on'      , ...
          'YMinorTick'  , 'on'      , ...
          'YGrid'       , 'on'      , ...
          'XColor'      , [0 0 0], ...
          'YColor'      , [0 0 0], ...
          'LineWidth'   , 1         );

        if strcmp(location,'X')
            % get limits and update limit text boxes
            xlimits = get(handles.axes_graphs,'XLim');
            ylimits = get(handles.axes_graphs,'YLim');
            set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
            set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
            set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
            set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
        end

        if strcmp(location,'n') % set new figure limits
            axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
        end

        if strcmp(location,'s')
            dir_str = strcat(handles.pathname,'plots', filesep); % in the plots dir
            if(~exist(dir_str, 'dir'))
                mkdir(dir_str);
            end
            file_str = strcat(dir_str,'plot_data_ppe.dat'); % create dat file with the data
            dlmwrite(file_str, ppe, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');

            file_str = strcat(dir_str,'plot_data_ppe_hist.dat'); % create dat file with the hist data
            dlmwrite(file_str, [xout' n'], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
        end
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  fun_update_listbox() % 

handles = guidata(gcbo);

    % handles.event_local_struct_info
    % 1 id
    % 2 integral
    % 3 FWHM_dwell
    % 4 number of local struct
    % 5 event start ind
    % 6 event end ind
    % 7 region integral
    % 8 max peak FWHM
    % 9 group 1 type:
    %       -1 - bad
    %       0 - unsorted
    %       1 - folded
    %       2 - clogged
    %       3 - unfolded
    %       
    % 10 group 2 type
    %       0 - unsorted
    %       1 - bare (no spikes)
    %       2 - spike A
    %       3 - spike B
    %       4 - spike C
    %       5 - spike D
    %       6 - spike E
    %       7 - spike F
    %       8 - spike G
    %       9 - spike H
    %       10 - spike I
    % 11 number of local struct from region analysis
    % 12 mean first extra points
    % 13 mean last extra points
    % 14 selected 1 -yes 0 - no

    % event types:
% 1 - All Events
% 2 - Current Selection
% 3 - Unsorted
% 4 - Folded
% 5 - Clogged
% 6 - Unfolded - All
% 7 - Unfolded - Unsorted
% 8 - Unfolded - Bare
% 9 - Unfolded - Spike A
% 10 Unfolded - Spike B
% 11 Unfolded - Spike C
% 12 Unfolded - Spike D
% 13 Unfolded - Spike E
% 14 Unfolded - Spike F
% 15 Unfolded - Spike G
% 16 Unfolded - Spike H
% 17 Unfolded - Spike I
% 18 Bad

switch get(handles.pop_select_group,'Value') % what to show
    case 1 % All Events
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info;
        
    case 2 % 2 - Current Selection
        
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,14)==1),:);
        
    case 3 % 3 - Unsorted
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==0),:);
        
    case 4 % 4 - Folded
        display('hello')
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==1),:);
        
        
    case 5 % 5 - Clogged
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==2),:);
        
        
    case 6 % 6 - Unfolded - All
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3),:);
        
        
    case 7 % 7 - Unfolded - Unsorted
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==0),:);
        
        
    case 8 % 8 - Unfolded - Bare
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==1),:);
        
        
    case 9 % 9 - Unfolded - Spike A
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==2),:);
        
        
    case 10 % 10 Unfolded - Spike B
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==3),:);
        
        
    case 11 % 11 Unfolded - Spike C
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==4),:);
        
        
    case 12 % 12 Unfolded - Spike D
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==5),:);
        
    case 13 % 13 Unfolded - Spike E
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==6),:);
        
        
    case 14 % 14 Unfolded - Spike F
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==7),:);
        
        
    case 15 % 15 Unfolded - Spike G
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==8),:);
        
        
    case 16 % 16 Unfolded - Spike H
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==9),:);
        
        
    case 17 % 17 Unfolded - Spike I
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==10),:);
        
        
    case 18 % 18 Bad
        
        handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==-1),:);
    
    
end

switch get(handles.pop_sort_events,'Value') % how to sort
    case 1 % ID
        
        handles.event_local_struct_info_listbox = sortrows(handles.event_local_struct_info_listbox,1);
        
    case 2 % # of Peaks
        
        handles.event_local_struct_info_listbox = sortrows(handles.event_local_struct_info_listbox,4);
        
    case 3 % Integral
        
        handles.event_local_struct_info_listbox = sortrows(handles.event_local_struct_info_listbox,2);        
        
    case 4 % Dwell time
        
        handles.event_local_struct_info_listbox = sortrows(handles.event_local_struct_info_listbox,3);
        
    otherwise
end

last_event_shown = size(handles.event_local_struct_info_listbox,1);
% assignin('base', 'last_event_shown', last_event_shown)
if get(handles.listbox_eventsbox,'Value') > last_event_shown
    set(handles.listbox_eventsbox,'Value',last_event_shown)
end

set(handles.listbox_eventsbox,'String',num2str(handles.event_local_struct_info_listbox))

total_events = size(handles.event_local_struct_info,1);
set(handles.edit_total_events,'String',num2str(total_events))
    
total_peaks = size(handles.peak_info,1);
set(handles.edit_total_peaks,'String',num2str(total_peaks))

listed_events = size(handles.event_local_struct_info_listbox,1);
set(handles.edit_selected_events,'String',num2str(listed_events))
    
% listed_peaks = size(find(handles.peak_info(:,1)==handles.event_local_struct_info_listbox(:,1),1);
% set(handles.edit_total_peaks,'String',num2str(listed_peaks))

% fun_plot_event('X', handles) % plot the event that was clicked

guidata(gcbo, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_conductance_Callback(~, ~, ~)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_save_event_Callback(~, ~, handles)


        handles.eventselected = get(handles.listbox_eventsbox,'Value'); % find event selected
        eventid = handles.event_local_struct_info_listbox(handles.eventselected,1);
        newFig = figure;
        fun_plot_event('s',handles) % data file saved inside this function
        saveas(newFig, strcat(handles.pathname, filesep,'events', filesep,'event_',num2str(eventid),'.fig'), 'fig'); % save as .fig
%         saveas(newFig, strcat(handles.pathname, filesep,'events', filesep,'event_',num2str(eventid),'.png'), 'png');
        set(newFig, 'PaperPositionMode', 'auto');
        saveas(newFig, strcat(handles.pathname, filesep,'events', filesep,'event_',num2str(eventid),'.png'), 'png');
%         print(newFig,'-depsc2',strcat(handles.pathname, filesep,'events', filesep,'event_',num2str(eventid),'.eps')) 

        if exist('export_fig') == 2
            set(gcf, 'Color', 'w');
            export_fig(newFig, strcat(handles.pathname, filesep,'events', filesep,'event_',num2str(eventid),'_expfig.eps'), '-cmyk');
            export_fig(newFig, strcat(handles.pathname, filesep,'events', filesep,'event_',num2str(eventid),'_expfig.png'), '-cmyk');
        end
        close(newFig)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_show_extrapoints_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_select_group_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_select_group_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_region_threshold_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_region_threshold_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_region_selection_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_region_selection_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_region_param_a_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_region_param_a_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_region_analyze_Callback(hObject, eventdata, handles)


set(handles.button_open,'BackgroundColor',[1 0 0])
indpos = 1; % pointer to current position, initialize
handles.region_analysis_peak_info = [];

extrapoints = handles.readextrapoints; % get number of extra points to each side of the event

newtrace = 0;
newstart = 0;

num_events_make_new_trace = size(handles.event_info,1);

% pre allocation
newstart = zeros(num_events_make_new_trace);
minmax_within_events = zeros(num_events_make_new_trace,10);
newend = zeros(num_events_make_new_trace);

paramA = str2double(get(handles.edit_region_param_a,'String'));
baselineDNA = str2double(get(handles.edit_baseline,'String'));

bottom_level = paramA*baselineDNA;

for i=1:num_events_make_new_trace % for number of events
    
    eventpoints = [(handles.event_info(i,8)) (handles.event_info(i,9))]; % get event start stop
    
    eventlength = handles.event_info(i,9) - handles.event_info(i,8) + 2*extrapoints + 1; % calc length in points
    eventtrace = handles.rawtrace_vect(eventpoints(1)-extrapoints:eventpoints(2)+extrapoints) - handles.event_info(i,5); % get the event trace, and remove baseline
%     if get(handles.check_flip,'Value')
%         eventtrace = eventtrace.*-1;
%     end
    mean_first_extra = abs(mean(eventtrace(1:extrapoints)));
    handles.event_local_struct_info(i,12) = mean_first_extra;
    mean_last_extra = abs(mean(eventtrace(end-extrapoints+1:end)));
    handles.event_local_struct_info(i,13) = mean_last_extra;
    mean_all_extra = abs(mean([mean_first_extra mean_last_extra]));
    
    eventtrace = eventtrace + baselineDNA; % get the event trace, and remove  first DNA level
    
    eventtrace(1:extrapoints) = zeros([extrapoints,1]);
    eventtrace(end-extrapoints+1:end) = zeros([extrapoints,1]);
%     if get(handles.check_flip,'Value')
%     	eventtrace = eventtrace.*-1;
%     end
    newstart(i) = indpos + extrapoints; % calculate the new start position
    
    eventtrace_inside = handles.rawtrace_vect(eventpoints(1):eventpoints(2)) - handles.event_info(i,5) + baselineDNA;
    timetemp = transpose(0:(handles.timestep*1000):handles.timestep*(eventpoints(2)-eventpoints(1))*1000); % build event time vector
    
%     if get(handles.check_flip,'Value')
%     	eventtrace_inside = eventtrace_inside.*-1;
%     end
    
    for j=1:eventlength - 2*extrapoints
        point_value = eventtrace_inside(j);
        if point_value > 0
            eventtrace_inside(j) = 0;
        elseif point_value < -1*bottom_level
            eventtrace_inside(j) = -1*bottom_level;
        end
    end
    
    eventtrace(extrapoints+1:end-extrapoints) = eventtrace_inside;
    
    event_int = -1.*trapz(timetemp,eventtrace_inside); % in nA*ms
    event_int_norm = event_int/(handles.event_local_struct_info(i,3)*bottom_level); %dimensionless
    handles.event_local_struct_info(i,7) = event_int_norm; % put region integral norm in matrix
    
    newtrace(indpos:indpos+eventlength-1) = eventtrace; % add this trace to new trace
    
    newend(i) = newstart(i) + handles.event_info(i,9) - handles.event_info(i,8); % calculate the new end position
    indpos = indpos + eventlength; % calc new start index value
end

handles.region_analysis_trace = newtrace;
guidata(gcbo, handles);
fun_find_local_struct(-1, 0, handles)
handles=guidata(gcbo); % get handles


%%%%%%%%%%%%%%%%%%%%%%%%
min_region_integral = min(handles.event_local_struct_info(:,7));
max_region_integral = max(handles.event_local_struct_info(:,7));
% set(handles.edit_classify_region_int_min,'String',num2str(min_region_integral))
set(handles.edit_classify_region_int_min,'String','0')
set(handles.edit_classify_region_int_max,'String',num2str(max_region_integral))
%%%%%%%%%%%%%%%%%%%%%%%%
% set(handles.edit_classify_FWHM_min,'String',num2str(min(handles.event_local_struct_info(:,8))))
set(handles.edit_classify_FWHM_min,'String','0')
set(handles.edit_classify_FWHM_max,'String',num2str(max(handles.event_local_struct_info(:,8))))
%%%%%%%%%%%%%%%%%%%%%%%%
min_region_numpeak = min(handles.event_local_struct_info(:,11));
max_region_numpeak = max(handles.event_local_struct_info(:,11));
set(handles.edit_classify_numpeak_min,'String',num2str(min_region_numpeak))
set(handles.edit_classify_numpeak_max,'String',num2str(max_region_numpeak))
%%%%%%%%%%%%%%%%%%%%%%%%
min_region_first_extra = min(handles.event_local_struct_info(:,12));
max_region_first_extra = max(handles.event_local_struct_info(:,12));
min_region_last_extra = min(handles.event_local_struct_info(:,13));
max_region_last_extra = max(handles.event_local_struct_info(:,13));
% set(handles.edit_first_extra_min_mean,'String',num2str(min_region_first_extra))
set(handles.edit_first_extra_min_mean,'String','0')
set(handles.edit_first_extra_max_mean,'String',num2str(max_region_first_extra))
% set(handles.edit_last_extra_min_mean,'String',num2str(min_region_last_extra))
set(handles.edit_last_extra_min_mean,'String','0')
set(handles.edit_last_extra_max_mean,'String',num2str(max_region_last_extra))
%%%%%%%%%%%%%%%%%%%%%%%%

set(handles.edit_region_peaks_found,'String',num2str(size(handles.region_analysis_peak_info,1)))

guidata(gcbo, handles);
fun_update_listbox()
set(handles.radio_region_int,'Value',1);
set(handles.radio_region_max_FWHM,'Value',0);
set(handles.radio_region_num_peaks,'Value',0);
set(handles.radio_region_max_mean_extra,'Value',0);
fun_plot_region_int(handles,'X') % plot
set(handles.button_open,'BackgroundColor',[0.925 0.914 0.847])
guidata(hObject, handles); % update handles



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_region_int_Callback(hObject, eventdata, handles)
% hObject    handle to radio_region_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_region_int
set(handles.radio_region_int,'Value',1);
set(handles.radio_region_max_FWHM,'Value',0);
set(handles.radio_region_num_peaks,'Value',0);
set(handles.radio_region_max_mean_extra,'Value',0);
set(handles.radio_region_max_mean_extra_last,'Value',0);
fun_plot_region_int(handles,'X') % plot
guidata(hObject, handles); % update handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_region_int(handles, location)
    colormap('default')
    if strcmp(location,'X') % plot here
        cla(handles.axes_graphs,'reset')
        axes(handles.axes_graphs)
    end
    if strcmp(location,'n') % plot in new figure
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        figure()
    end

    bins = str2num(get(handles.edit_classify_region_int_bins,'String'));
    
    intdata = handles.event_local_struct_info(:,7);
    

    [n,xout] = hist(intdata,bins); % plot it
    hist(intdata,bins)
    binwidth = xout(2)-xout(1);
    set(handles.edit_2_bin_width,'String',num2str(binwidth))

    hXLabel = xlabel('Normalized Integral between I_1 and I_2');
    hYLabel = ylabel('Counts');
    
    % style
    set( gca                       , ...
        'FontName'   , 'Arial' );
    set([hXLabel, hYLabel], ...
        'FontName'   , 'Arial');
    set([hXLabel, hYLabel]  , ...
        'FontSize'   , 10          );
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'on'      , ...
      'YMinorTick'  , 'on'      , ...
      'YGrid'       , 'on'      , ...
      'XColor'      , [0 0 0], ...
      'YColor'      , [0 0 0], ...
      'LineWidth'   , 1         );
  
    if strcmp(location,'X')
        % get limits and update limit text boxes
        xlimits = get(handles.axes_graphs,'XLim');
        ylimits = get(handles.axes_graphs,'YLim');
        set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
        set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
        set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
        set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
    end

    if strcmp(location,'n')% use same limits
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
    end

    if strcmp(location,'s')
        dir_str = strcat(handles.pathname,'plots', filesep); % in the plots directory
        if(~exist(dir_str, 'dir'))
            mkdir(dir_str);
        end
        file_str = strcat(dir_str,'plot_data_region_int.dat'); % save dat file with dwell data
        dlmwrite(file_str, intdata, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
        
        file_str = strcat(dir_str,'plot_data_region_int_hist.dat'); % create dat file with the hist data
        dlmwrite(file_str, [xout' n'], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    end
    if strcmp(location,'n')% normalized cumulative plot
        if exist('export_fig') == 2
            set(gcf, 'Color', 'w');
            export_fig(gcf, strcat(handles.pathname,'plots', filesep,'plot_norm_hist_region.eps'), '-cmyk');
            export_fig(gcf, strcat(handles.pathname,'plots', filesep,'plot_norm_hist_region.png'), '-cmyk');
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this is for the cumulative figure

    
    n_sum = cumsum(n);
    n_sum_norm = n_sum./n_sum(end);
    
%     frac_folded_events = 0.64;
%     frac_unfolded_events = 1 - frac_folded_events;
%     index_vert = find(n_sum_norm > frac_unfolded_events);
%     vert_intersect_point = xout(index_vert(1));

    figure()
    bar(xout,n_sum_norm)
    hold on
%     line([0 1],[frac_unfolded_events frac_unfolded_events],'Color','red')
%     line([vert_intersect_point vert_intersect_point],[0 1],'Color','red')
    hXLabel = xlabel('Normalized Integral between I_1 and I_2');
    hYLabel = ylabel('Normalized Cumulative Sum');
    
    if strcmp(location,'n')% normalized cumulative plot
        if exist('export_fig') == 2
            set(gcf, 'Color', 'w');
            export_fig(gcf, strcat(handles.pathname,'plots', filesep,'plot_norm_cumsum_region.eps'), '-cmyk');
            export_fig(gcf, strcat(handles.pathname,'plots', filesep,'plot_norm_cumsum_region.png'), '-cmyk');
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_region_max_FWHM_Callback(hObject, eventdata, handles)

set(handles.radio_region_int,'Value',0);
set(handles.radio_region_max_FWHM,'Value',1);
set(handles.radio_region_num_peaks,'Value',0);
set(handles.radio_region_max_mean_extra,'Value',0);
set(handles.radio_region_max_mean_extra_last,'Value',0);
fun_plot_region_max_FWHM(handles,'X') % plot
guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_region_max_FWHM(handles, location)
    colormap('default')
    if strcmp(location,'X') % plot here
        cla(handles.axes_graphs,'reset')
        axes(handles.axes_graphs)
    end
    if strcmp(location,'n') % plot in new figure
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        figure()
    end

    bins = str2num(get(handles.edit_classify_FWHM_bins,'String'));
    
    % use max region integral value to select spikes?
    use_max_region_int = get(handles.check_use_max_region_int_in_FWHM,'Value');
    
    if use_max_region_int
        
        max_region_integral = str2double(get(handles.edit_classify_region_int_max,'String'));
        
        eventid_region_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,7)<max_region_integral),1);
        index_region_peaks_to_analyze = find(ismember(handles.region_analysis_peak_info(:,1),eventid_region_index_selected));
        
        FWHMdata = handles.region_analysis_peak_info(index_region_peaks_to_analyze,5);
        
    else
        % spike dwell data
        FWHMdata = handles.region_analysis_peak_info(:,5);
    
    end
    % 1st DNA level crossing
    % FWHMdata = handles.region_analysis_peak_info_base_cross_dwell;

    [n,xout] = hist(FWHMdata,bins); % plot it
    hist(FWHMdata,bins)
    binwidth = xout(2)-xout(1);
    set(handles.edit_2_bin_width,'String',num2str(binwidth))

    hXLabel = xlabel('FWHM (ms)');
    hYLabel = ylabel('Counts');
    
    % style
    set( gca                       , ...
        'FontName'   , 'Arial' );
    set([hXLabel, hYLabel], ...
        'FontName'   , 'Arial');
    set([hXLabel, hYLabel]  , ...
        'FontSize'   , 10          );
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'on'      , ...
      'YMinorTick'  , 'on'      , ...
      'YGrid'       , 'on'      , ...
      'XColor'      , [0 0 0], ...
      'YColor'      , [0 0 0], ...
      'LineWidth'   , 1         );
  
    if strcmp(location,'X')
        % get limits and update limit text boxes
        xlimits = get(handles.axes_graphs,'XLim');
        ylimits = get(handles.axes_graphs,'YLim');
        set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
        set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
        set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
        set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
    end

    if strcmp(location,'n')% use same limits
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
    end

    if strcmp(location,'s')
        dir_str = strcat(handles.pathname,'plots', filesep); % in the plots directory
        if(~exist(dir_str, 'dir'))
            mkdir(dir_str);
        end
        file_str = strcat(dir_str,'plot_data_region_FWHM.dat'); % save dat file with dwell data
        dlmwrite(file_str, FWHMdata, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
        
        file_str = strcat(dir_str,'plot_data_region_FWHM_hist.dat'); % create dat file with the hist data
        dlmwrite(file_str, [xout' n'], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_region_num_peaks_Callback(hObject, eventdata, handles)

set(handles.radio_region_int,'Value',0);
set(handles.radio_region_max_FWHM,'Value',0);
set(handles.radio_region_num_peaks,'Value',1);
set(handles.radio_region_max_mean_extra,'Value',0);
set(handles.radio_region_max_mean_extra_last,'Value',0);
fun_plot_region_num_peaks(handles,'X') % plot
guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_region_num_peaks(handles, location)
    colormap('default')
    if strcmp(location,'X') % plot here
        cla(handles.axes_graphs,'reset')
        axes(handles.axes_graphs)
    end
    if strcmp(location,'n') % plot in new figure
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        figure()
    end

    bins = str2num(get(handles.edit_classify_numpeak_bins,'String'));
    
    numpeakdata = handles.event_local_struct_info(:,11);
    

    [n,xout] = hist(numpeakdata,bins); % plot it
    hist(numpeakdata,bins)
    binwidth = xout(2)-xout(1);
    set(handles.edit_2_bin_width,'String',num2str(binwidth))

    hXLabel = xlabel('Integral (nA*ms)');
    hYLabel = ylabel('Counts');
    
    % style
    set( gca                       , ...
        'FontName'   , 'Arial' );
    set([hXLabel, hYLabel], ...
        'FontName'   , 'Arial');
    set([hXLabel, hYLabel]  , ...
        'FontSize'   , 10          );
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'on'      , ...
      'YMinorTick'  , 'on'      , ...
      'YGrid'       , 'on'      , ...
      'XColor'      , [0 0 0], ...
      'YColor'      , [0 0 0], ...
      'LineWidth'   , 1         );
  
    if strcmp(location,'X')
        % get limits and update limit text boxes
        xlimits = get(handles.axes_graphs,'XLim');
        ylimits = get(handles.axes_graphs,'YLim');
        set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
        set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
        set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
        set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
    end

    if strcmp(location,'n')% use same limits
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
    end

    if strcmp(location,'s')
        dir_str = strcat(handles.pathname,'plots', filesep); % in the plots directory
        if(~exist(dir_str, 'dir'))
            mkdir(dir_str);
        end
        file_str = strcat(dir_str,'plot_data_region_numpeak.dat'); % save dat file with dwell data
        dlmwrite(file_str, numpeakdata, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
        
        file_str = strcat(dir_str,'plot_data_region_numpeak_hist.dat'); % create dat file with the hist data
        dlmwrite(file_str, [xout' n'], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_region_max_mean_extra_Callback(hObject, eventdata, handles)

set(handles.radio_region_int,'Value',0);
set(handles.radio_region_max_FWHM,'Value',0);
set(handles.radio_region_num_peaks,'Value',0);
set(handles.radio_region_max_mean_extra,'Value',1);
set(handles.radio_region_max_mean_extra_last,'Value',0);
fun_plot_region_max_mean_extra_first(handles,'X') % plot
guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_region_max_mean_extra_first(handles, location)
    colormap('default')
    if strcmp(location,'X') % plot here
        cla(handles.axes_graphs,'reset')
        axes(handles.axes_graphs)
    end
    if strcmp(location,'n') % plot in new figure
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        figure()
    end

    bins = str2num(get(handles.edit_region_bins_mean_extra,'String'));
    


    [~,xout] = hist(handles.event_local_struct_info(:,12),bins); % plot it
    
%     subplot(2,1,1,'Parent',handles.axes_graphs); 
    
    hist(handles.event_local_struct_info(:,12),bins)
    binwidth = xout(2)-xout(1);
    set(handles.edit_2_bin_width,'String',num2str(binwidth))

    hXLabel = xlabel('Mean Current Value (nA)');
    hYLabel = ylabel('Counts');
    
    % style
    set( gca                       , ...
        'FontName'   , 'Arial' );
    set([hXLabel, hYLabel], ...
        'FontName'   , 'Arial');
    set([hXLabel, hYLabel]  , ...
        'FontSize'   , 10          );
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'on'      , ...
      'YMinorTick'  , 'on'      , ...
      'YGrid'       , 'on'      , ...
      'XColor'      , [0 0 0], ...
      'YColor'      , [0 0 0], ...
      'LineWidth'   , 1         );
  
%     subplot(2,1,2,'Parent',handles.axes_graphs); 
% 
%     [n2,xout2] = hist(last_extra_data,bins); % plot it
%    
%     
%     hist(last_extra_data,bins)
% 
% 
%     hXLabel = xlabel('Mean Current Value (nA)');
%     hYLabel = ylabel('Counts');
    
    if strcmp(location,'X')
        % get limits and update limit text boxes
        xlimits = get(handles.axes_graphs,'XLim');
        ylimits = get(handles.axes_graphs,'YLim');
        set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
        set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
        set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
        set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
    end

    if strcmp(location,'n')% use same limits
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
    end

    if strcmp(location,'s')
        dir_str = strcat(handles.pathname,'plots', filesep); % in the plots directory
        if(~exist(dir_str, 'dir'))
            mkdir(dir_str);
        end
%         file_str = strcat(dir_str,'plot_data_region_int.dat'); % save dat file with dwell data
%         dlmwrite(file_str, intdata, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
%         
%         file_str = strcat(dir_str,'plot_data_region_int_hist.dat'); % create dat file with the hist data
%         dlmwrite(file_str, [xout' n'], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_region_max_mean_extra_last(handles, location)
    colormap('default')
    if strcmp(location,'X') % plot here
        cla(handles.axes_graphs,'reset')
        axes(handles.axes_graphs)
    end
    if strcmp(location,'n') % plot in new figure
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        figure()
    end

    bins = str2num(get(handles.edit_region_bins_mean_extra,'String'));
    

%     last_extra_data = [handles.event_local_struct_info(:,12) handles.event_local_struct_info(:,13)];

    [n,xout] = hist(handles.event_local_struct_info(:,13),bins); % plot it
    
%     subplot(2,1,1,'Parent',handles.axes_graphs); 
    
    hist(handles.event_local_struct_info(:,13),bins)
    binwidth = xout(2)-xout(1);
    set(handles.edit_2_bin_width,'String',num2str(binwidth))

    hXLabel = xlabel('Mean Current Value (nA)');
    hYLabel = ylabel('Counts');
    
    % style
    set( gca                       , ...
        'FontName'   , 'Arial' );
    set([hXLabel, hYLabel], ...
        'FontName'   , 'Arial');
    set([hXLabel, hYLabel]  , ...
        'FontSize'   , 10          );
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'on'      , ...
      'YMinorTick'  , 'on'      , ...
      'YGrid'       , 'on'      , ...
      'XColor'      , [0 0 0], ...
      'YColor'      , [0 0 0], ...
      'LineWidth'   , 1         );
  
%     subplot(2,1,2,'Parent',handles.axes_graphs); 
% 
%     [n2,xout2] = hist(last_extra_data,bins); % plot it
%    
%     
%     hist(last_extra_data,bins)
% 
% 
%     hXLabel = xlabel('Mean Current Value (nA)');
%     hYLabel = ylabel('Counts');
    
    if strcmp(location,'X')
        % get limits and update limit text boxes
        xlimits = get(handles.axes_graphs,'XLim');
        ylimits = get(handles.axes_graphs,'YLim');
        set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
        set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
        set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
        set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
    end

    if strcmp(location,'n')% use same limits
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
    end

    if strcmp(location,'s')
        dir_str = strcat(handles.pathname,'plots', filesep); % in the plots directory
        if(~exist(dir_str, 'dir'))
            mkdir(dir_str);
        end
%         file_str = strcat(dir_str,'plot_data_region_int.dat'); % save dat file with dwell data
%         dlmwrite(file_str, intdata, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
%         
%         file_str = strcat(dir_str,'plot_data_region_int_hist.dat'); % create dat file with the hist data
%         dlmwrite(file_str, [xout' n'], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_region_int_bins_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_region_int_bins_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_region_int_min_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_region_int_min_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_region_int_max_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_region_int_max_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_FWHM_bins_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_FWHM_bins_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_FWHM_min_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_FWHM_min_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_FWHM_max_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_FWHM_max_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_numpeak_bins_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_numpeak_bins_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_numpeak_min_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_numpeak_min_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_numpeak_max_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_classify_numpeak_max_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_select_Callback(hObject, eventdata, handles)
% hObject    handle to button_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.event_local_struct_info(:,14) = zeros(size(handles.event_local_struct_info,1),1); % reset selected

min_region_integral = str2double(get(handles.edit_classify_region_int_min,'String'));
max_region_integral = str2double(get(handles.edit_classify_region_int_max,'String'));

min_region_FWHM = str2double(get(handles.edit_classify_FWHM_min,'String'));
max_region_FWHM = str2double(get(handles.edit_classify_FWHM_max,'String'));

min_region_numpeak = str2double(get(handles.edit_classify_numpeak_min,'String'));
max_region_numpeak = str2double(get(handles.edit_classify_numpeak_max,'String'));

region_first_extra_min_mean = str2double(get(handles.edit_first_extra_min_mean,'String'));
region_first_extra_max_mean = str2double(get(handles.edit_first_extra_max_mean,'String'));

region_last_extra_min_mean = str2double(get(handles.edit_last_extra_min_mean,'String'));
region_last_extra_max_mean = str2double(get(handles.edit_last_extra_max_mean,'String'));

% int 7, max FWHM 8, num peaks region 11, max extra 12 & 13
index_selected = find(handles.event_local_struct_info(:,7) >= min_region_integral & ...
                      handles.event_local_struct_info(:,7) <= max_region_integral & ...
                      handles.event_local_struct_info(:,8) >= min_region_FWHM & ...
                      handles.event_local_struct_info(:,8) <= max_region_FWHM & ...
                      handles.event_local_struct_info(:,11) >= min_region_numpeak & ...
                      handles.event_local_struct_info(:,11) <= max_region_numpeak & ...
                      handles.event_local_struct_info(:,12) >= region_first_extra_min_mean & ...
                      handles.event_local_struct_info(:,12) <= region_first_extra_max_mean & ...
                      handles.event_local_struct_info(:,13) >= region_last_extra_min_mean & ...
                      handles.event_local_struct_info(:,13) <= region_last_extra_max_mean & ...
                      handles.event_local_struct_info(:,9) == 0);

handles.event_local_struct_info(index_selected,14) = ones(size(index_selected));

set(handles.edit_num_events_selected_classify,'String',num2str(size(index_selected,1)))

set(handles.pop_select_group,'Value',2) % show selected

guidata(hObject, handles); % update handles

fun_update_listbox()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_assign_group_Callback(hObject, eventdata, handles)
% hObject    handle to button_assign_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


index_selected = find(handles.event_local_struct_info(:,14)==1);

% Unsorted
% Folded
% Clogged
% Unfolded - Unsorted
% Unfolded - Bare
% Unfolded - Spike A
% Unfolded - Spike B
% Unfolded - Spike C
% Unfolded - Spike D
% Unfolded - Spike E
% Unfolded - Spike F
% Unfolded - Spike G
% Unfolded - Spike H
% Unfolded - Spike I
% Bad

switch get(handles.pop_assign_group,'Value') % what to show
    case 1 % Unsorted
        handles.event_local_struct_info(index_selected,9) = 0;
        
    case 2 % Folded
        handles.event_local_struct_info(index_selected,9) = 1;
        
    case 3 % Clogged
        handles.event_local_struct_info(index_selected,9) = 2;
        
    case 4 % Unfolded - Unsorted
        handles.event_local_struct_info(index_selected,9) = 3;
        handles.event_local_struct_info(index_selected,10) = 0;
        
    case 5 % Unfolded - Bare
        handles.event_local_struct_info(index_selected,9) = 3;
         handles.event_local_struct_info(index_selected,10) = 1;
        
    case 6 % Unfolded - Spike A
        handles.event_local_struct_info(index_selected,9) = 3;
         handles.event_local_struct_info(index_selected,10) = 2;
        
    case 7 % Unfolded - Spike B
        handles.event_local_struct_info(index_selected,9) = 3;
         handles.event_local_struct_info(index_selected,10) = 3;
        
    case 8 % Unfolded - Spike C
        handles.event_local_struct_info(index_selected,9) = 3;
         handles.event_local_struct_info(index_selected,10) = 4;
        
    case 9 % Unfolded - Spike D
        handles.event_local_struct_info(index_selected,9) = 3;
         handles.event_local_struct_info(index_selected,10) = 5;
        
    case 10 % Unfolded - Spike E
        handles.event_local_struct_info(index_selected,9) = 3;
         handles.event_local_struct_info(index_selected,10) = 6;
        
    case 11 % Unfolded - Spike F
        handles.event_local_struct_info(index_selected,9) = 3;
         handles.event_local_struct_info(index_selected,10) = 7;
        
    case 12 % Unfolded - Spike G
        handles.event_local_struct_info(index_selected,9) = 3;
         handles.event_local_struct_info(index_selected,10) = 8;
        
    case 13 % Unfolded - Spike H
        handles.event_local_struct_info(index_selected,9) = 3;
         handles.event_local_struct_info(index_selected,10) = 9;
        
    case 14 % Unfolded - Spike I
        handles.event_local_struct_info(index_selected,9) = 3;
         handles.event_local_struct_info(index_selected,10) = 10;
        
    case 15 % Bad
        handles.event_local_struct_info(index_selected,9) = -1;


end

guidata(gcbo, handles); % save handles

fun_update_listbox()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_assign_group_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_assign_group_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_assign_event_Callback(hObject, eventdata, handles)
% hObject    handle to button_assign_event (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.check_use_ids,'Value') == 0
    handles.eventselected = get(handles.listbox_eventsbox,'Value'); % find event selected

    eventid = handles.event_local_struct_info_listbox(handles.eventselected,1);
    

    % Unsorted
    % Folded
    % Clogged
    % Unfolded - Unsorted
    % Unfolded - Bare
    % Unfolded - Spike A
    % Unfolded - Spike B
    % Unfolded - Spike C
    % Unfolded - Spike D
    % Unfolded - Spike E
    % Unfolded - Spike F
    % Unfolded - Spike G
    % Unfolded - Spike H
    % Unfolded - Spike I
    % Bad

    switch get(handles.pop_assign_event,'Value') % what to show
        case 1 % Unsorted
            handles.event_local_struct_info(eventid,9) = 0;

        case 2 % Folded
            handles.event_local_struct_info(eventid,9) = 1;

        case 3 % Clogged
            handles.event_local_struct_info(eventid,9) = 2;

        case 4 % Unfolded - Unsorted
            handles.event_local_struct_info(eventid,9) = 3;
            handles.event_local_struct_info(eventid,10) = 0;

        case 5 % Unfolded - Bare
            handles.event_local_struct_info(eventid,9) = 3;
             handles.event_local_struct_info(eventid,10) = 1;

        case 6 % Unfolded - Spike A
            handles.event_local_struct_info(eventid,9) = 3;
             handles.event_local_struct_info(eventid,10) = 2;

        case 7 % Unfolded - Spike B
            handles.event_local_struct_info(eventid,9) = 3;
             handles.event_local_struct_info(eventid,10) = 3;

        case 8 % Unfolded - Spike C
            handles.event_local_struct_info(eventid,9) = 3;
             handles.event_local_struct_info(eventid,10) = 4;

        case 9 % Unfolded - Spike D
            handles.event_local_struct_info(eventid,9) = 3;
             handles.event_local_struct_info(eventid,10) = 5;

        case 10 % Unfolded - Spike E
            handles.event_local_struct_info(eventid,9) = 3;
             handles.event_local_struct_info(eventid,10) = 6;

        case 11 % Unfolded - Spike F
            handles.event_local_struct_info(eventid,9) = 3;
             handles.event_local_struct_info(eventid,10) = 7;

        case 12 % Unfolded - Spike G
            handles.event_local_struct_info(eventid,9) = 3;
             handles.event_local_struct_info(eventid,10) = 8;

        case 13 % Unfolded - Spike H
            handles.event_local_struct_info(eventid,9) = 3;
             handles.event_local_struct_info(eventid,10) = 9;

        case 14 % Unfolded - Spike I
            handles.event_local_struct_info(eventid,9) = 3;
             handles.event_local_struct_info(eventid,10) = 10;

        case 15 % Bad
            handles.event_local_struct_info(eventid,9) = -1;


    end
end
if get(handles.check_use_ids,'Value') == 1
    
    id_array = (str2num(get(handles.edit_resort_ids,'String')));
    number_of_events_to_resort = size(id_array,1);
    
    if number_of_events_to_resort > 0
        for ijk = 1:number_of_events_to_resort
            
            switch get(handles.pop_assign_event,'Value') % what to show
                case 1 % Unsorted
                    handles.event_local_struct_info(id_array(ijk),9) = 0;

                case 2 % Folded
                    handles.event_local_struct_info(id_array(ijk),9) = 1;

                case 3 % Clogged
                    handles.event_local_struct_info(id_array(ijk),9) = 2;

                case 4 % Unfolded - Unsorted
                    handles.event_local_struct_info(id_array(ijk),9) = 3;
                    handles.event_local_struct_info(id_array(ijk),10) = 0;

                case 5 % Unfolded - Bare
                    handles.event_local_struct_info(id_array(ijk),9) = 3;
                     handles.event_local_struct_info(id_array(ijk),10) = 1;

                case 6 % Unfolded - Spike A
                    handles.event_local_struct_info(id_array(ijk),9) = 3;
                     handles.event_local_struct_info(id_array(ijk),10) = 2;

                case 7 % Unfolded - Spike B
                    handles.event_local_struct_info(id_array(ijk),9) = 3;
                     handles.event_local_struct_info(id_array(ijk),10) = 3;

                case 8 % Unfolded - Spike C
                    handles.event_local_struct_info(id_array(ijk),9) = 3;
                     handles.event_local_struct_info(id_array(ijk),10) = 4;

                case 9 % Unfolded - Spike D
                    handles.event_local_struct_info(id_array(ijk),9) = 3;
                     handles.event_local_struct_info(id_array(ijk),10) = 5;

                case 10 % Unfolded - Spike E
                    handles.event_local_struct_info(id_array(ijk),9) = 3;
                     handles.event_local_struct_info(id_array(ijk),10) = 6;

                case 11 % Unfolded - Spike F
                    handles.event_local_struct_info(id_array(ijk),9) = 3;
                     handles.event_local_struct_info(id_array(ijk),10) = 7;

                case 12 % Unfolded - Spike G
                    handles.event_local_struct_info(id_array(ijk),9) = 3;
                     handles.event_local_struct_info(id_array(ijk),10) = 8;

                case 13 % Unfolded - Spike H
                    handles.event_local_struct_info(id_array(ijk),9) = 3;
                     handles.event_local_struct_info(id_array(ijk),10) = 9;

                case 14 % Unfolded - Spike I
                    handles.event_local_struct_info(id_array(ijk),9) = 3;
                     handles.event_local_struct_info(id_array(ijk),10) = 10;

                case 15 % Bad
                    handles.event_local_struct_info(id_array(ijk),9) = -1;


            end
            
        end
        
    end    
    
end

guidata(gcbo, handles); % save handles
fun_update_listbox()
handles=guidata(gcbo); % get handles
fun_plot_event('X', handles) % plot the next event

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_gen_selected_spikes_list()

handles = guidata(gcbo);

if size(handles.peak_info,1) > 0
    switch get(handles.pop_analyse_peaks_from,'Value') % what to show
        case 1 % All Events

            handles.eventid_index_selected = handles.event_local_struct_info(:,1);
            handles.index_peaks_to_analyze = handles.peak_info(:,2);

        case 2 % 2 - Current Selection

            handles.eventid_index_selected = find(handles.event_local_struct_info(:,14)==1);
%             handles.index_peaks_to_analyze = handles.event_local_struct_info(find(handles.peak_info(:,1)==handles.eventid_index_selected),:);
            

        case 3 % 3 - Unsorted

            handles.eventid_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==0),1);
            
            

        case 4 % 4 - Folded

            handles.eventid_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==1),1);


        case 5 % 5 - Clogged

            handles.eventid_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==2),1);


        case 6 % 6 - Unfolded - All

            handles.eventid_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3),1);


        case 7 % 7 - Unfolded - Unsorted

            handles.eventid_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==0),1);


        case 8 % 8 - Unfolded - Bare

            handles.eventid_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==1),1);


        case 9 % 9 - Unfolded - Spike A

            handles.eventid_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==2),1);


        case 10 % 10 Unfolded - Spike B

            handles.eventid_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==3),1);


        case 11 % 11 Unfolded - Spike C

            handles.eventid_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==4),1);


        case 12 % 12 Unfolded - Spike D

            handles.eventid_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==5),1);

        case 13 % 13 Unfolded - Spike E

            handles.eventid_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==6),1);


        case 14 % 14 Unfolded - Spike F

            handles.eventid_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==7),1);


        case 15 % 15 Unfolded - Spike G

            handles.eventid_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==8),1);


        case 16 % 16 Unfolded - Spike H

            handles.event_local_struct_info_listbox = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==9),1);


        case 17 % 17 Unfolded - Spike I

            handles.eventid_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==10),1);


        case 18 % 18 Bad

            handles.eventid_index_selected = handles.event_local_struct_info(find(handles.event_local_struct_info(:,9)==-1),1);


    end
%     assignin('base', 'peak_info', handles.peak_info)
%     assignin('base', 'eventid_index_selected', handles.eventid_index_selected)
    handles.index_peaks_to_analyze = find(ismember(handles.peak_info(:,1),handles.eventid_index_selected));
 
%     assignin('base', 'index_peaks_to_analyze', handles.index_peaks_to_analyze)
end

if isempty(handles.index_peaks_to_analyze)
    handles.index_peaks_to_analyze = -1;
end
if isempty(handles.eventid_index_selected)
    handles.eventid_index_selected = -1;
end

guidata(gcbo, handles); % save handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_assign_event_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_assign_event_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_show_unfiltered_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function buttom_dump_var_Callback(hObject, eventdata, handles)


% this can be useful to rename variables for specific voltages
use_voltage_in_names = 0;

baselineDNA = str2double(get(handles.edit_baseline,'String'));

if use_voltage_in_names
    
    assignin('base', strcat('event_info_',num2str(handles.voltage),'mV'), handles.event_info)
    assignin('base', strcat('peak_info_',num2str(handles.voltage),'mV'), handles.peak_info)
    assignin('base', strcat('region_analysis_peak_info_',num2str(handles.voltage),'mV'), handles.region_analysis_peak_info)
    if isfield(handles,'event_more_info') % only do this if it exists
        assignin('base', strcat('event_more_info_',num2str(handles.voltage),'mV'), handles.event_more_info)
    end
    assignin('base', strcat('readextrapoints_',num2str(handles.voltage),'mV'), handles.readextrapoints)
    assignin('base', strcat('rawtrace_vect_',num2str(handles.voltage),'mV'), handles.rawtrace_vect)
    assignin('base', strcat('event_local_struct_info_',num2str(handles.voltage),'mV'), handles.event_local_struct_info)
    assignin('base', strcat('region_analysis_peak_info_',num2str(handles.voltage),'mV'), handles.region_analysis_peak_info)
    assignin('base', strcat('DNA_baseline_',num2str(handles.voltage),'mV'), baselineDNA)
    
else
    
    assignin('base', 'event_info', handles.event_info)
    assignin('base', 'peak_info', handles.peak_info)
    assignin('base', 'region_analysis_peak_info', handles.region_analysis_peak_info)
    if isfield(handles,'event_more_info') % only do this if it exists
        assignin('base', 'event_more_info', handles.event_more_info)
    end
    assignin('base', 'readextrapoints', handles.readextrapoints)
    assignin('base', 'rawtrace_vect', handles.rawtrace_vect)
    assignin('base', 'event_local_struct_info', handles.event_local_struct_info)
    assignin('base', 'region_analysis_peak_info', handles.region_analysis_peak_info)
    assignin('base', 'DNA_baseline', baselineDNA)
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_X_events_calc_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_X_events_calc_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_first_X_mean_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_first_X_mean_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_first_X_STD_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_first_X_STD_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_last_X_mean_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_last_X_mean_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_last_X_STD_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_last_X_STD_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_selected_events_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_selected_events_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_selected_peaks_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_selected_peaks_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_analyse_peaks_from_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_analyse_peaks_from_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pos_abs_bins_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pos_abs_bins_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_show_region_analysis_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_region_peaks_found_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_region_peaks_found_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_num_events_selected_classify_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_num_events_selected_classify_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_region_bins_mean_extra_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_region_bins_mean_extra_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_first_extra_max_mean_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_first_extra_max_mean_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_last_extra_max_mean_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_last_extra_max_mean_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_show_thresholds_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pop_dist_Callback(hObject, eventdata, handles)

fun_plot_pop_dist(handles,'X')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_pop_dist(handles, location)

        if strcmp(location,'X') % plot here
            cla(handles.axes_graphs,'reset')
            axes(handles.axes_graphs)
        end
        if strcmp(location,'n') % plot in new figure
            ylimits = get(handles.axes_graphs,'YLim');
            xlimits = get(handles.axes_graphs,'XLim');
            figure()
        end

% -1 - bad
    %       0 - unsorted
    %       1 - folded
    %       2 - clogged
    %       3 - unfolded
    %       
    % 10 group 2 type
    %       0 - unsorted
    %       1 - bare (no spikes)
    %       2 - spike A
    %       3 - spike B
    %       4 - spike C
    %       5 - spike D
    %       6 - spike E
    %       7 - spike F
    %       8 - spike G
    %       9 - spike H
    %       10 - spike I
    
num_bad = size(find(handles.event_local_struct_info(:,9)==-1),1);
num_unsorted = size(find(handles.event_local_struct_info(:,9)==0),1);
num_folded = size(find(handles.event_local_struct_info(:,9)==1),1);
num_clogged = size(find(handles.event_local_struct_info(:,9)==2),1);
num_unfolded_unsorted = size(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==0),1);
num_unfolded_bare = size(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==1),1);
num_unfolded_sA = size(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==2),1);
num_unfolded_sB = size(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==3),1);
num_unfolded_sC = size(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==4),1);
num_unfolded_sD = size(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==5),1);
num_unfolded_sE = size(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==6),1);
num_unfolded_sF = size(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==7),1);
num_unfolded_sG = size(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==8),1);
num_unfolded_sH = size(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==9),1);
num_unfolded_sI = size(find(handles.event_local_struct_info(:,9)==3 & handles.event_local_struct_info(:,10)==10),1);
    
x = [num_unsorted num_folded num_clogged num_unfolded_unsorted num_unfolded_bare num_unfolded_sA num_unfolded_sB ...
     num_unfolded_sC num_unfolded_sD num_unfolded_sE num_unfolded_sF num_unfolded_sG num_unfolded_sH num_unfolded_sI];

num_unsorted_label = strcat('Unsorted (',num2str(num_unsorted),') [Bad (',num2str(num_bad),')]');
num_folded_label = strcat('Folded (',num2str(num_folded),')');
num_clogged_label = strcat('Clogged (',num2str(num_clogged),')');
num_unfolded_unsorted_label = strcat('Unfolded - Unsorted (',num2str(num_unfolded_unsorted),')');
num_unfolded_bare_label = strcat('Unfolded - Bare (',num2str(num_unfolded_bare),')');
num_unfolded_sA_label = strcat('Unfolded - Spike A (',num2str(num_unfolded_sA),')');
num_unfolded_sB_label = strcat('Unfolded - Spike B (',num2str(num_unfolded_sB),')');
num_unfolded_sC_label = strcat('Unfolded - Spike C (',num2str(num_unfolded_sC),')');
num_unfolded_sD_label = strcat('Unfolded - Spike D (',num2str(num_unfolded_sD),')');
num_unfolded_sE_label = strcat('Unfolded - Spike E (',num2str(num_unfolded_sE),')');
num_unfolded_sF_label = strcat('Unfolded - Spike F (',num2str(num_unfolded_sF),')');
num_unfolded_sG_label = strcat('Unfolded - Spike G (',num2str(num_unfolded_sG),')');
num_unfolded_sH_label = strcat('Unfolded - Spike H (',num2str(num_unfolded_sH),')');
num_unfolded_sI_label = strcat('Unfolded - Spike I (',num2str(num_unfolded_sI),')');

% labels = {'Unsorted () [Bad ()]','Folded','Clogged','Unfolded - Unsorted ()', ...
%     'Unfolded - Bare ()','Unfolded - Spike A ()','Unfolded - Spike B ()','Unfolded - Spike C ()',...
%     'Unfolded - Spike D ()','Unfolded - Spike E ()','Unfolded - Spike F ()','Unfolded - Spike G ()',...
%     'Unfolded - Spike H ()','Unfolded - Spike I ()'};

labels = {num_unsorted_label, num_folded_label, num_clogged_label, num_unfolded_unsorted_label, num_unfolded_bare_label, num_unfolded_sA_label, num_unfolded_sB_label, ...
     num_unfolded_sC_label, num_unfolded_sD_label, num_unfolded_sE_label, num_unfolded_sF_label, num_unfolded_sG_label, num_unfolded_sH_label, num_unfolded_sI_label};
 
% cla(handles.axes_graphs,'reset')
% axes(handles.axes_graphs);
h = pie(x,labels);

str=get(h(2:2:end), 'String');
assignin('base', 'str', str);
AX=legend(h(1:2:end), cat(1, str),'Location','EastOutside');
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',8);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_region_max_mean_extra_last_Callback(hObject, eventdata, handles)

set(handles.radio_region_int,'Value',0);
set(handles.radio_region_max_FWHM,'Value',0);
set(handles.radio_region_num_peaks,'Value',0);
set(handles.radio_region_max_mean_extra,'Value',0);
set(handles.radio_region_max_mean_extra_last,'Value',1);
fun_plot_region_max_mean_extra_last(handles,'X') % plot
guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_region_bins_mean_extra_last_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_region_bins_mean_extra_last_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_first_extra_min_mean_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_first_extra_min_mean_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_last_extra_min_mean_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_last_extra_min_mean_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_reset_region_int_Callback(hObject, eventdata, handles)

min_region_integral = min(handles.event_local_struct_info(:,7));
max_region_integral = max(handles.event_local_struct_info(:,7));

set(handles.edit_classify_region_int_min,'String',num2str(min_region_integral))
set(handles.edit_classify_region_int_max,'String',num2str(max_region_integral))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_reset_region_FWHM_Callback(hObject, eventdata, handles)

set(handles.edit_classify_FWHM_min,'String',num2str(min(handles.event_local_struct_info(:,8))))
set(handles.edit_classify_FWHM_max,'String',num2str(max(handles.event_local_struct_info(:,8))))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_reset_region_numpeaks_Callback(hObject, eventdata, handles)

min_region_numpeak = min(handles.event_local_struct_info(:,11));
max_region_numpeak = max(handles.event_local_struct_info(:,11));

set(handles.edit_classify_numpeak_min,'String',num2str(min_region_numpeak))
set(handles.edit_classify_numpeak_max,'String',num2str(max_region_numpeak))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_reset_region_first_ex_mean_Callback(hObject, eventdata, handles)

min_region_first_extra = min(handles.event_local_struct_info(:,12));
max_region_first_extra = max(handles.event_local_struct_info(:,12));
set(handles.edit_first_extra_min_mean,'String',num2str(min_region_first_extra))
set(handles.edit_first_extra_max_mean,'String',num2str(max_region_first_extra))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_reset_region_last_ex_mean_Callback(hObject, eventdata, handles)

min_region_last_extra = min(handles.event_local_struct_info(:,13));
max_region_last_extra = max(handles.event_local_struct_info(:,13));

set(handles.edit_last_extra_min_mean,'String',num2str(min_region_last_extra))
set(handles.edit_last_extra_max_mean,'String',num2str(max_region_last_extra))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_pos_abs_Callback(hObject, eventdata, handles)

set(handles.radio_scatter_spike,'Value',0);
set(handles.radio_dwell,'Value',0);
set(handles.radio_amp,'Value',0);
set(handles.radio_pos,'Value',0);
set(handles.radio_pos_abs,'Value',1);
set(handles.radio_current,'Value',0);
set(handles.radio_peaks_per_event,'Value',0);
fun_plot_pos_abs(handles, 'X')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_pos_abs(handles, location)


fun_gen_selected_spikes_list()
    handles = guidata(gcbo);
    
    if handles.index_peaks_to_analyze ~= -1

        colormap('default')
        if strcmp(location,'X') % plot here
            cla(handles.axes_graphs,'reset')
            axes(handles.axes_graphs)
        end
        if strcmp(location,'n') % plot in new figure
            ylimits = get(handles.axes_graphs,'YLim');
            xlimits = get(handles.axes_graphs,'XLim');
            figure()
        end

        bins = str2num(get(handles.edit_pos_abs_bins,'String'));

        posdata = handles.peak_info(handles.index_peaks_to_analyze,4);


        [n,xout] = hist(posdata,bins); % plot it
        hist(posdata,bins)
        binwidth = xout(2)-xout(1);
        set(handles.edit_2_bin_width,'String',num2str(binwidth))

        hXLabel = xlabel('Temporal Position (ms)');
        hYLabel = ylabel('Counts');

        % style
        set( gca                       , ...
            'FontName'   , 'Arial' );
        set([hXLabel, hYLabel], ...
            'FontName'   , 'Arial');
        set([hXLabel, hYLabel]  , ...
            'FontSize'   , 10          );
        set(gca, ...
          'Box'         , 'on'     , ...
          'TickDir'     , 'in'     , ...
          'TickLength'  , [.02 .02] , ...
          'XMinorTick'  , 'on'      , ...
          'YMinorTick'  , 'on'      , ...
          'YGrid'       , 'on'      , ...
          'XColor'      , [0 0 0], ...
          'YColor'      , [0 0 0], ...
          'LineWidth'   , 1         );

        if strcmp(location,'X')
            % get limits and update limit text boxes
            xlimits = get(handles.axes_graphs,'XLim');
            ylimits = get(handles.axes_graphs,'YLim');
            set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
            set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
            set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
            set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
        end

        if strcmp(location,'n')% use same limits
            axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
        end

        if strcmp(location,'s')
            dir_str = strcat(handles.pathname,'plots', filesep); % in the plots directory
            if(~exist(dir_str, 'dir'))
                mkdir(dir_str);
            end
            file_str = strcat(dir_str,'plot_data_posabs.dat'); % save dat file with dwell data
            dlmwrite(file_str, posdata, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');

            file_str = strcat(dir_str,'plot_data_posabs_hist.dat'); % create dat file with the hist data
            dlmwrite(file_str, [xout' n'], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
        end
    end
    
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_calc_first_last_X_events_ppe_mean_STD_Callback(hObject, eventdata, handles)

numevents = size(handles.event_local_struct_info(:,1),1);

first_last_X_events = str2double(get(handles.edit_X_events_calc,'String'));

if numevents >= 2*first_last_X_events && size(handles.peak_info(:,1),1) > 0
    
    first_X_ppe = handles.event_local_struct_info(1:first_last_X_events,4);
    
    last_X_ppe = handles.event_local_struct_info(end-first_last_X_events+1:end,4);
    
    set(handles.edit_first_X_mean,'String',num2str(mean(first_X_ppe)))
    set(handles.edit_first_X_STD,'String',num2str(std(first_X_ppe)))
    set(handles.edit_last_X_mean,'String',num2str(mean(last_X_ppe)))
    set(handles.edit_last_X_STD,'String',num2str(std(last_X_ppe)))
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_reset_listbox_Callback(hObject, eventdata, handles)

set(handles.listbox_eventsbox,'Value',1)
guidata(gcbo, handles); % update handles



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pk_stat_total_events_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pk_stat_total_events_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pk_stat_events_with_spike_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pk_stat_events_with_spike_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pk_stat_percent_of_ends_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pk_stat_percent_of_ends_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pk_stat_percent_events_start_spike_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pk_stat_percent_events_start_spike_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pk_stat_percent_events_end_spike_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pk_stat_percent_events_end_spike_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_update_peak_stats

    fun_gen_selected_spikes_list()
    % handles.index_peaks_to_analyze
    % handles.eventid_index_selected
    
    handles=guidata(gcbo); % get handles
    
    
    percent_from_ends = str2num(get(handles.edit_pk_stat_percent_of_ends,'String'))/100;
    number_of_events = size(handles.eventid_index_selected,1);
    ppe = handles.event_local_struct_info(handles.eventid_index_selected,4);
    max_ppe = max(ppe);
    set(handles.edit_ppe_bins,'String',num2str(max_ppe))
    total_spikes = sum(ppe);
    index_events_with_spikes = find(ppe>=1);
    events_with_spikes = size(index_events_with_spikes,1);
    events_with_spikes_start = 0;
    events_with_spikes_end = 0;
    
    number_of_spikes_start = 0;
    number_of_spikes_middle = 0;
    number_of_spikes_end = 0;
    
    
    
%     assignin('base', 'peak_info', handles.peak_info)
%     assignin('base', 'ppe', ppe)
%     assignin('base', 'eventid_index_selected', handles.eventid_index_selected)
    
    if number_of_events > 1
        if ~isempty(handles.index_peaks_to_analyze)
            
            for i=1:number_of_events
                
                index_current_peaks = find(handles.peak_info(:,1)==handles.eventid_index_selected(i));
                
%                 assignin('base', 'index_current_peaks', index_current_peaks)
                
                if ~isempty(index_current_peaks)
                    
                    pos_values = handles.peak_info(index_current_peaks,10);
                    
%                     disp(strcat(num2str(i),' here ',num2str(size(pos_values,2))))
                    
                    start_flag = 0;
                    end_flag = 0;
                    
                    for j=1:size(pos_values,1)
                        if pos_values(j) <= percent_from_ends
                            start_flag = 1;
                            number_of_spikes_start = number_of_spikes_start + 1;
                        end
                        if pos_values(j) >= (1-percent_from_ends)
                            end_flag = 1;
                            number_of_spikes_end = number_of_spikes_end + 1;
                        end
                        if pos_values(j) >= (0.5-percent_from_ends/2) && pos_values(j) <= (0.5+percent_from_ends/2)
                            number_of_spikes_middle = number_of_spikes_middle + 1;
                        end
                    end
                    if start_flag == 1
                        events_with_spikes_start = events_with_spikes_start+1;
                    end
                    if end_flag == 1
                        events_with_spikes_end = events_with_spikes_end+1;
                    end
                end
            end
            
        end
    end
    
    pcent_spikes_start = number_of_spikes_start/total_spikes;
    pcent_spikes_middle = number_of_spikes_middle/total_spikes;
    pcent_spikes_end = number_of_spikes_end/total_spikes;
    set(handles.edit_spikes_at_start_pcent,'String',num2str(pcent_spikes_start))
    set(handles.edit_spikes_in_middle_pcent,'String',num2str(pcent_spikes_middle))
    set(handles.edit_spikes_at_end_pcent,'String',num2str(pcent_spikes_end))
    
    set(handles.edit_pk_stat_total_events,'String',num2str(number_of_events))
    set(handles.edit_pk_stat_events_with_spike,'String',num2str(events_with_spikes))
    set(handles.edit_selected_peaks,'String',num2str(total_spikes))
    set(handles.edit_pk_stat_percent_events_start_spike,'String',num2str(events_with_spikes_start/number_of_events))
    set(handles.edit_pk_stat_percent_events_end_spike,'String',num2str(events_with_spikes_end/number_of_events))
    
    guidata(gcbo, handles); % save handles
    
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_fit_gaussian_Callback(hObject, eventdata, handles)

    fun_gen_selected_spikes_list()
    handles = guidata(gcbo);
    
%     disp('here 1')
    
    if handles.index_peaks_to_analyze ~= -1
        
%         disp('here 2')

        HandleMainGUI=getappdata(0,'HandleMainGUI'); % handle used to pass data to new GUI

        dwell_vector = handles.peak_info(handles.index_peaks_to_analyze,5);
        
        setappdata(HandleMainGUI,'dwellvector',dwell_vector); % pass peak dwell data to sub GUI
        
        fun_make_new_trace(handles) % reformat trace so that baselines are zero and cutoff extrapoints if necessary
        handles=guidata(gcbo); % update local handles
        
        setappdata(HandleMainGUI,'newtrace',handles.newtrace); % pass it down

        setappdata(HandleMainGUI,'dwellbins',str2num(get(handles.edit_dwell_bins,'String'))); % pass number of dwell bins
        setappdata(HandleMainGUI,'cbbins',str2num(get(handles.edit_amp_bins,'String'))); % pass number of amplitude bins
        setappdata(HandleMainGUI,'currentbins',str2num(get(handles.edit_current_bins,'String'))); % pass number of current bins

        setappdata(HandleMainGUI,'amplitude',handles.peak_info(handles.index_peaks_to_analyze,8)); % pass it down

        setappdata(HandleMainGUI,'check_conductance',0); % use current

        setappdata(HandleMainGUI,'pathname',handles.pathname);
        GUI_gaussian % open gaussian fitting GUI
    end

    guidata(hObject, handles); % update handles
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_time_Callback(hObject, eventdata, handles)

    fun_plot_time(handles, 'n')
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_time(handles, location)
    fun_gen_selected_spikes_list()
    handles = guidata(gcbo);
    
    if handles.index_peaks_to_analyze ~= -1

        colormap('default')
        if strcmp(location,'X') % plot in GUI
            cla(handles.axes_graphs,'reset')
            axes(handles.axes_graphs)
        end

        if strcmp(location,'n') % plot in new figure
%             ylimits = get(handles.axes_graphs,'YLim');
%             xlimits = get(handles.axes_graphs,'XLim');
            figure()
        end

        bins = str2num(get(handles.edit_ppe_bins,'String'));

        index_of_events_with_peaks = find(handles.event_local_struct_info(:,4)>0);
        index_of_selected_events_with_peaks = index_of_events_with_peaks(ismember(index_of_events_with_peaks,handles.eventid_index_selected));
        
%         assignin('base', 'eventid_index_selected', handles.eventid_index_selected)
%         
%         assignin('base', 'index_of_events_with_peaks', index_of_events_with_peaks)
%         assignin('base', 'index_of_selected_events_with_peaks', index_of_selected_events_with_peaks)
        
        index_of_events_with_no_peaks = find(handles.event_local_struct_info(:,4)==0);
        index_of_selected_events_with_no_peaks = index_of_events_with_no_peaks(ismember(index_of_events_with_no_peaks,handles.eventid_index_selected));
        
%         assignin('base', 'index_of_events_with_no_peaks', index_of_selected_events_with_no_peaks)
%         assignin('base', 'index_of_selected_events_with_no_peaks', index_of_selected_events_with_no_peaks)
%         
%         assignin('base', 'event_more_info', handles.event_more_info)
        
        time_vect_peaks = handles.event_more_info(index_of_selected_events_with_peaks,4);
        time_vect_no_peaks = handles.event_more_info(index_of_selected_events_with_no_peaks,4);

        [n_peak,xout_peak] = hist(time_vect_peaks,bins); % plot it
        [n_no_peak,xout_no_peak] = hist(time_vect_no_peaks,xout_peak); % plot it
        
        percentage_of_events_with_spike = (n_peak./(n_peak+n_no_peak));
        
        mean_percentage = sum(n_peak)/(sum(n_peak)+sum(n_no_peak));
        
        std_spikes = std(percentage_of_events_with_spike);
        
        display(strcat('time mean: ',num2str(mean_percentage)))
        display(strcat('time STD: ',num2str(std_spikes)))
        
        total_events = (sum(n_peak)+sum(n_no_peak));
        
%         assignin('base', 'n_peak', n_peak)
%         assignin('base', 'n_no_peak', n_no_peak)
%         assignin('base', 'xout_peak', xout_peak)
        bar(xout_peak,percentage_of_events_with_spike)
        axis([0 max(xout_peak) 0 1])
        
%         binwidth = xout(2)-xout(1);
%         set(handles.edit_2_bin_width,'String',num2str(binwidth))

    %     if get(handles.check_conductance,'Value')
    %         hXLabel = xlabel('Conductance (nS)');
    %     else
            hXLabel = xlabel('Time (s)');
    %     end

        hYLabel = ylabel('Fraction of Events with Spike');

        % style
        set( gca                       , ...
            'FontName'   , 'Arial' );
        set([hXLabel, hYLabel], ...
            'FontName'   , 'Arial');
        set([hXLabel, hYLabel]  , ...
            'FontSize'   , 10          );

        set(gca, ...
          'Box'         , 'on'     , ...
          'TickDir'     , 'in'     , ...
          'TickLength'  , [.02 .02] , ...
          'XMinorTick'  , 'on'      , ...
          'YMinorTick'  , 'on'      , ...
          'YGrid'       , 'on'      , ...
          'XColor'      , [0 0 0], ...
          'YColor'      , [0 0 0], ...
          'LineWidth'   , 1         );

        if strcmp(location,'X')
            % get limits and update limit text boxes
            xlimits = get(handles.axes_graphs,'XLim');
            ylimits = get(handles.axes_graphs,'YLim');
            set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
            set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
            set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
            set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
        end

        if strcmp(location,'n') % set new figure limits
%             axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
        end

        if strcmp(location,'s') % s
            dir_str = strcat(handles.pathname,'plots', filesep); % in the plots dir
            if(~exist(dir_str, 'dir'))
                mkdir(dir_str);
            end
            file_str = strcat(dir_str,'plot_data_time_spike.dat'); % create dat file with the data
            dlmwrite(file_str, time_vect_peaks, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
            
            file_str = strcat(dir_str,'plot_data_time_no_spike.dat'); % create dat file with the data
            dlmwrite(file_str, time_vect_no_peaks, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');

            file_str = strcat(dir_str,'plot_data_time_hist.dat'); % create dat file with the hist data
            dlmwrite(file_str, [xout_peak' percentage_of_events_with_spike'], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
        end
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_spikes_at_start_pcent_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_spikes_at_start_pcent_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_spikes_in_middle_pcent_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_spikes_in_middle_pcent_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_spikes_at_end_pcent_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_spikes_at_end_pcent_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_scatter_bins_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_scatter_bins_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_scatter_spike_Callback(hObject, eventdata, handles)

set(handles.radio_scatter_spike,'Value',1);
set(handles.radio_dwell,'Value',0);
set(handles.radio_amp,'Value',0);
set(handles.radio_pos,'Value',0);
set(handles.radio_pos_abs,'Value',0);
set(handles.radio_current,'Value',0);
set(handles.radio_peaks_per_event,'Value',0);
fun_plot_scatter_spike(handles,'X') % plot current hist
guidata(hObject, handles); % update handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_scatter_spike(handles, location)


    fun_gen_selected_spikes_list()
    handles = guidata(gcbo);
    
    if handles.index_peaks_to_analyze ~= -1

        colormap('default')
        if strcmp(location,'X') % plot here
            cla(handles.axes_graphs,'reset')
            axes(handles.axes_graphs)
        end
        if strcmp(location,'n') % plot in new figure
            ylimits = get(handles.axes_graphs,'YLim');
            xlimits = get(handles.axes_graphs,'XLim');
            figure()
        end

        
        amplitude = handles.peak_info(handles.index_peaks_to_analyze,8);
        
        if get(handles.check_conductance,'Value')
            amplitude = amplitude./(handles.voltage/1000);
        end
        
        use_pos_instead_of_dwell = 0;
        
        if use_pos_instead_of_dwell
            dwelldata = handles.peak_info(handles.index_peaks_to_analyze,10); % position data
        else
            dwelldata = handles.peak_info(handles.index_peaks_to_analyze,5); % dwell data
        end

%         assignin('base', 'amplitude', amplitude)
%         assignin('base', 'dwelldata', dwelldata)
        
        scatter_hist_on = 1;
        
        if scatter_hist_on && strcmp(location,'n') % scatter hist plot, called by mass export button press
            
            pp3 = scatterhist(dwelldata,amplitude,[str2num(get(handles.edit_dwell_bins,'String')) str2num(get(handles.edit_amp_bins,'String'))]);
        
        else % otherwise use normal scatter plot

            pp3 = scatter(dwelldata,amplitude,'or');

        end
        
        if use_pos_instead_of_dwell
            hXLabel = xlabel('Normalized Position');
            
        else
            hXLabel = xlabel('Dwell Time (ms)');
            
        end
        
        
        if get(handles.check_conductance,'Value')
            hYLabel = ylabel('Conductance (nS)');
        else
            hYLabel = ylabel('Current Blockade (nA)');
        end

        % style
        set( gca                       , ...
            'FontName'   , 'Arial' );
        set([hXLabel, hYLabel], ...
            'FontName'   , 'Arial');
        set([hXLabel, hYLabel]  , ...
            'FontSize'   , 10          );
        set(gca, ...
          'Box'         , 'on'     , ...
          'TickDir'     , 'in'     , ...
          'TickLength'  , [.02 .02] , ...
          'XMinorTick'  , 'on'      , ...
          'YMinorTick'  , 'on'      , ...
          'YGrid'       , 'on'      , ...
          'XColor'      , [0 0 0], ...
          'YColor'      , [0 0 0], ...
          'LineWidth'   , 1         );

        if strcmp(location,'X')
            % get limits and update limit text boxes
            xlimits = get(handles.axes_graphs,'XLim');
            ylimits = get(handles.axes_graphs,'YLim');
            set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
            set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
            set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
            set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
        end

        if strcmp(location,'n')% use same limits
            axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
        end
    end

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_region_plot_Callback(hObject, eventdata, handles)

if get(handles.radio_region_int,'Value')
    fun_plot_region_int(handles, 'n')
end
if get(handles.radio_region_max_FWHM,'Value')
    fun_plot_region_max_FWHM(handles, 'n')
end
if get(handles.radio_region_num_peaks,'Value')
    fun_plot_region_num_peaks(handles, 'n')
end
if get(handles.radio_region_max_mean_extra,'Value')
    fun_plot_region_max_mean_extra_first(handles, 'n')
end
if get(handles.radio_region_max_mean_extra_last,'Value')
    fun_plot_region_max_mean_extra_last(handles, 'n')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_use_max_region_int_in_FWHM_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_basline_in_FWHM_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_save_all_events_Callback(hObject, eventdata, handles)

total_events = size(handles.event_local_struct_info_listbox,1);

dir_str = strcat(handles.pathname, filesep,'all_events', filesep); % in the all_events subdir
if(~exist(dir_str, 'dir'))
    mkdir(dir_str);
end

for i=1:total_events
        set(handles.listbox_eventsbox,'Value',i); % set event
        eventselected  = i;
        eventid = handles.event_local_struct_info_listbox(eventselected,1);
        newFig = figure;
        fun_plot_event('z',handles) % only plot
        set(newFig, 'PaperPositionMode', 'auto');
        saveas(newFig, strcat(handles.pathname, filesep,'all_events', filesep,'event_',num2str(eventid),'.png'), 'png');
        close(newFig)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_peakdwelldata_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_peakdwelldata_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_resort_partof_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_resort_partof_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_resort_eventhasnumpeaks_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_resort_eventhasnumpeaks_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_resort_assignto_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_resort_assignto_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_resort_events_Callback(hObject, eventdata, handles)

% All Events
% Current Selection
% Unsorted
% Folded
% Clogged
% Unfolded - Unsorted
% Unfolded - Bare
% Unfolded - Spike A
% Unfolded - Spike B
% Unfolded - Spike C
% Unfolded - Spike D
% Unfolded - Spike E
% Unfolded - Spike F
% Unfolded - Spike G
% Unfolded - Spike H
% Unfolded - Spike I
% Bad

    % handles.event_local_struct_info
    % 9 group 1 type:
    %       -1 - bad
    %       0 - unsorted
    %       1 - folded
    %       2 - clogged
    %       3 - unfolded
    %       
    % 10 group 2 type
    %       0 - unsorted
    %       1 - bare (no spikes)
    %       2 - spike A
    %       3 - spike B
    %       4 - spike C
    %       5 - spike D
    %       6 - spike E
    %       7 - spike F
    %       8 - spike G
    %       9 - spike H
    %       10 - spike I

resort_group_from = get(handles.pop_resort_partof,'Value');
resort_number_peaks = str2num(get(handles.edit_resort_eventhasnumpeaks,'String'));
resort_group_to = get(handles.pop_resort_assignto,'Value');

switch resort_group_from % where the events are from 
    case 1 % All Events
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,4) ==  resort_number_peaks);
        
    case 2 % Current Selection
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,14) ==  1 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
        
    case 3 % Unsorted
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  0 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
        
    case 4 % Folded
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  1 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
        
    case 5 % Clogged
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  2 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
       
    case 6 % Unfolded - All
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  3 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
    
    case 7 % Unfolded - Unsorted
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  3 & handles.event_local_struct_info(:,10) ==  0 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
        
    case 8 % Unfolded - Bare
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  2 & handles.event_local_struct_info(:,10) ==  1 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
        
    case 9 % Unfolded - Spike A
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  2 & handles.event_local_struct_info(:,10) ==  2 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
        
    case 10 % Unfolded - Spike B
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  2 & handles.event_local_struct_info(:,10) ==  3 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
        
    case 11 % Unfolded - Spike C
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  2 & handles.event_local_struct_info(:,10) ==  4 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
        
    case 12 % Unfolded - Spike D
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  2 & handles.event_local_struct_info(:,10) ==  5 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
        
    case 13 % Unfolded - Spike E
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  2 & handles.event_local_struct_info(:,10) ==  6 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
        
    case 14 % Unfolded - Spike F
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  2 & handles.event_local_struct_info(:,10) ==  7 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
        
    case 15 % Unfolded - Spike G
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  2 & handles.event_local_struct_info(:,10) ==  8 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
        
    case 16 % Unfolded - Spike H
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  2 & handles.event_local_struct_info(:,10) ==  9 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
        
    case 17 % Unfolded - Spike I
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  2 & handles.event_local_struct_info(:,10) ==  10 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);
        
    case 18 % Bad
        
        index_of_events_to_resort = find(handles.event_local_struct_info(:,9) ==  -1 & handles.event_local_struct_info(:,4) ==  resort_number_peaks);

end

switch resort_group_to % where the events moving to 
    case 1 % Unsorted
        
        handles.event_local_struct_info(index_of_events_to_resort,9) = 0;
        
    case 2 % Folded
        
        handles.event_local_struct_info(index_of_events_to_resort,9) = 1;
        
    case 3 % Clogged
        
        handles.event_local_struct_info(index_of_events_to_resort,9) = 2;
        
    case 4 % Unfolded - Unsorted
        
        handles.event_local_struct_info(index_of_events_to_resort,9) = 3;
        handles.event_local_struct_info(index_of_events_to_resort,10) = 0;
        
    case 5 % Unfolded - Bare
        
        handles.event_local_struct_info(index_of_events_to_resort,9) = 3;
        handles.event_local_struct_info(index_of_events_to_resort,10) = 1;
        
    case 6 % Unfolded - Spike A
        
        handles.event_local_struct_info(index_of_events_to_resort,9) = 3;
        handles.event_local_struct_info(index_of_events_to_resort,10) = 2;
        
    case 7 % Unfolded - Spike B
        
        handles.event_local_struct_info(index_of_events_to_resort,9) = 3;
        handles.event_local_struct_info(index_of_events_to_resort,10) = 3;
        
    case 8 % Unfolded - Spike C
        
        handles.event_local_struct_info(index_of_events_to_resort,9) = 3;
        handles.event_local_struct_info(index_of_events_to_resort,10) = 4;
        
    case 9 % Unfolded - Spike D
        
        handles.event_local_struct_info(index_of_events_to_resort,9) = 3;
        handles.event_local_struct_info(index_of_events_to_resort,10) = 5;
        
    case 10 % Unfolded - Spike E
        
        handles.event_local_struct_info(index_of_events_to_resort,9) = 3;
        handles.event_local_struct_info(index_of_events_to_resort,10) = 6;
        
    case 11 % Unfolded - Spike F
        
        handles.event_local_struct_info(index_of_events_to_resort,9) = 3;
        handles.event_local_struct_info(index_of_events_to_resort,10) = 7;
        
    case 12 % Unfolded - Spike G
        
        handles.event_local_struct_info(index_of_events_to_resort,9) = 3;
        handles.event_local_struct_info(index_of_events_to_resort,10) = 8;
        
    case 13 % Unfolded - Spike H
        
        handles.event_local_struct_info(index_of_events_to_resort,9) = 3;
        handles.event_local_struct_info(index_of_events_to_resort,10) = 9;
        
    case 14 % Unfolded - Spike I
        
        handles.event_local_struct_info(index_of_events_to_resort,9) = 3;
        handles.event_local_struct_info(index_of_events_to_resort,10) = 10;
        
    case 15 % Bad
        
        handles.event_local_struct_info(index_of_events_to_resort,9) = -1;
        
end

guidata(gcbo, handles); % save handles
fun_update_listbox()



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_resort_ids_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_resort_ids_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_use_ids_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_fwhm_center_peak_pos_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_FWHM_level_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_FWHM_level_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_list_peak_center_time_Callback(hObject, eventdata, handles)
    % empty
