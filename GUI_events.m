function varargout = GUI_events(varargin)
% GUI_EVENTS M-file for GUI_events.fig
%      GUI_EVENTS, by itself, creates a new GUI_EVENTS or raises the existing
%      singleton*.
%
%      H = GUI_EVENTS returns the handle to a new GUI_EVENTS or the handle
%      to
%      the existing singleton*.
%
%      GUI_EVENTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_EVENTS.M with the given input arguments.
%
%      GUI_EVENTS('Property','Value',...) creates a new GUI_EVENTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_events_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_events_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help GUI_events

% Last Modified by GUIDE v2.5 23-Oct-2014 09:15:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_events_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_events_OutputFcn, ...
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
function GUI_events_OpeningFcn(hObject, ~, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to GUI_events (see VARARGIN)

    % Choose default command line output for GUI_events
    handles.output = hObject;

    % MY INIVariables
    handles.gui_events_file_columns = [1 2 3 4 5 6 7 8 9]; % legacy variable
    handles.extrapoints = 50; % default number of extrapoints
    handles.timestep = 0.000002; % default time step, updated using sample rate in parameters file

    % 1average	2maximum	3integral	4dwelltime	5baseline	6detectionlevel
    % 7eventtype	8start point	9end point

    % set to box mode on startup
    set(handles.radio_ellipse_mode,'Value',0);
    set(handles.radio_box_mode,'Value',1);

    setappdata(0,'HandleMainGUI',hObject);

    % Update handles structure
    guidata(hObject, handles); % update handles

    % UIWAIT makes GUI_events wait for user response (see UIRESUME)
    % uiwait(handles.GUI_events);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Outputs from this function are returned to the command line.
function varargout = GUI_events_OutputFcn(~, ~, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);

    % Get default command line output from handles structure
    varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text_file_box_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text_file_box_CreateFcn(hObject, eventdata, handles)

    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%          listbox_eventsbox_Callback            %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if you click on an event, plot it
function listbox_eventsbox_Callback(hObject, ~, handles)
    % Hints: contents = get(hObject,'String') returns listbox_eventsbox contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from listbox_eventsbox
    set(handles.text_status,'String','Event selected...plotting...')
    if isfield(handles,'rawtrace_vect') % if rawtrace is loaded
        fun_plot_event('X') % plot the event that was clicked
    end
    set(handles.text_status,'String','Event selected...plotting...done')
    guidata(hObject, handles); % update handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                   fun_plot_event               %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this will plot a single event trace in the bottom left plot of the GUI
%
% Inputs:
%         location - 'X' to plot here in GUI
%                    'n' to plot in a new figure
%                    's' to have the event data
function fun_plot_event(location) % plot events in the GUI
    handles=guidata(gcbo); % get handles
    handles.eventselected = get(handles.listbox_eventsbox,'Value'); % find event selected
    set(handles.text_event_selected,'String',num2str(handles.eventselected)) % update event selected number

    % where did this event occur, here we figure out the file, segment, and
    % time
    full_location_number = handles.event_more_info_temp(handles.eventselected,1);
    file_number = floor(full_location_number/1000000);
    segment_number = mod(full_location_number,1000000);
    event_time = handles.event_more_info_temp(handles.eventselected,3);
    event_location_full_string = strcat('File:',{' '},num2str(file_number),{' '},'Seg:',{' '},num2str(segment_number),{' '},'Time:',{' '},num2str(event_time));
    
    set(handles.edit_event_location,'String',event_location_full_string) % update event location
        
    % get start stop points
    eventpoints = [(handles.event_info_temp(handles.eventselected,8)) (handles.event_info_temp(handles.eventselected,9))]; 
    
    % number of extrapoints
    handles.extrapoints = str2num(get(handles.edit_extra_points,'String')); 
    
    if get(handles.check_use_extra_points,'Value') % show extra points
        
        if get(handles.check_remove_baseline,'Value')
            eventtrace = handles.rawtrace_vect(eventpoints(1)-handles.extrapoints:eventpoints(2)+handles.extrapoints) - handles.event_info_temp(handles.eventselected,5); % build event trace
        else
            eventtrace = handles.rawtrace_vect(eventpoints(1)-handles.extrapoints:eventpoints(2)+handles.extrapoints); % build event trace
        end
        
        % build event time vector
        timetemp = transpose(0:(handles.timestep*1000):handles.timestep*(eventpoints(2)-eventpoints(1)+2*handles.extrapoints)*1000); 
        
        % show unfiltered event
        if get(handles.check_overlay_unfiltered,'Value') && isfield(handles,'rawunfilteredtrace_vect') 
            eventtraceunfilt = handles.rawunfilteredtrace_vect(eventpoints(1)-handles.extrapoints:eventpoints(2)+handles.extrapoints); % unfiltered trace
        end
        
    else % don't show extra points
        
        if get(handles.check_remove_baseline,'Value')
            % build event trace
            eventtrace = handles.rawtrace_vect(eventpoints(1):eventpoints(2)) - handles.event_info_temp(handles.eventselected,5); 
        
        else
            % build event trace
            eventtrace = handles.rawtrace_vect(eventpoints(1):eventpoints(2));
            
        end
        % build event time vector
        timetemp = transpose(0:(handles.timestep*1000):handles.timestep*(eventpoints(2)-eventpoints(1))*1000); 
        if get(handles.check_overlay_unfiltered,'Value') && isfield(handles,'rawunfilteredtrace_vect') % show unfiltered event?
            % unfiltered trace
            eventtraceunfilt = handles.rawunfilteredtrace_vect(eventpoints(1):eventpoints(2)); 
        end
    end

    if get(handles.check_conductance,'Value') % convert to nS?
        % get voltage
        voltage = str2double(get(handles.edit_voltage,'String'))/1000; 
        if get(handles.check_remove_baseline,'Value')
            % divide by voltage
            eventtrace = eventtrace./voltage; 
        else
            % remove baseline and divide by voltage
            eventtrace = (eventtrace- handles.event_info_temp(handles.eventselected,5))./voltage; 
        end
        if get(handles.check_overlay_unfiltered,'Value') && isfield(handles,'rawunfilteredtrace_vect') % show unfiltered event?
            % remove baseline and divide by voltage
            eventtraceunfilt = (eventtraceunfilt-handles.event_info_temp(handles.eventselected,5))./voltage; 
        end
    end
    if strcmp(location, 'X') % plot it here?
        cla(handles.axes_eventplot,'reset')
        axes(handles.axes_eventplot)
    else % or plot in new window
        newFigure = figure(); 
    end
    if get(handles.check_overlay_unfiltered,'Value') && isfield(handles,'rawunfilteredtrace_vect') % show unfiltered event?
        % plot unfiltered in cyan
        plot(timetemp,eventtraceunfilt,'-c') 
        hold on
    end
    if get(handles.check_show_baseline,'Value') % show baseline?
        if get(handles.check_conductance,'Value') || get(handles.check_remove_baseline,'Value')
            % in the conductance case it is zero
            baseline = zeros(length(eventtrace)); 
        else
            % otherwise read baseline from event info
            baseline = ones(length(eventtrace))*handles.event_info_temp(handles.eventselected,5); 
        end
        % plot baseline
        plot(timetemp,baseline,'-g') 
    end
    hold on
    % plot event trace
    pp4 = plot(timetemp,eventtrace,'-b'); 


    ylimits = get(handles.axes_eventplot,'YLim');
    xlimits = get(handles.axes_eventplot,'XLim');
    % scale to fit in X
    axis([0 max(timetemp) ylimits(1) ylimits(2)]) 

    hXLabel = xlabel('Time (ms)');
    if get(handles.check_conductance,'Value') % set proper Y label
        hYLabel = ylabel('Conductance (nS)');
    else
        hYLabel = ylabel('Current (nA)');
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
    set(pp4                     , ...
      'Color'           , [0 0 1], ...
      'LineWidth'       , 2           );
    
    % save data and plot
    if strcmp(location, 's')
        % in the plots dir
        dir_str = strcat(handles.pathname,'plots', filesep); 
        if(~exist(dir_str, 'dir')) % make dir if it doesn't exist
            mkdir(dir_str);
        end
        % create dat file with the data
        file_str = strcat(dir_str,'event_',num2str(handles.eventselected),'_filtered.dat'); 
        dlmwrite(file_str, eventtrace, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
        if get(handles.check_overlay_unfiltered,'Value') && isfield(handles,'rawunfilteredtrace_vect') % save unfiltered data
            % create dat file with the data
            file_str2 = strcat(dir_str,'event_',num2str(handles.eventselected),'_unfiltered.dat'); 
            dlmwrite(file_str2, eventtraceunfilt, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
        end
        
        %save figures in various formats
        saveas(newFigure, strcat(handles.pathname,'plots', filesep,num2str(handles.eventselected),'.fig'), 'fig'); % save as .fig
        saveas(newFigure, strcat(handles.pathname,'plots', filesep,num2str(handles.eventselected),'_large.png'), 'png');
        set(newFigure, 'PaperPositionMode', 'auto');
        saveas(newFigure, strcat(handles.pathname,'plots', filesep,num2str(handles.eventselected),'_small.png'), 'png');
        print(newFigure,'-depsc2',strcat(handles.pathname,'plots', filesep,num2str(handles.eventselected),'.eps')) 
        
        if exist('export_fig') == 2 % if export_fig exists save nice figures
            set(gcf, 'Color', 'w');
            export_fig(newFigure, strcat(handles.pathname,'plots', filesep,num2str(handles.eventselected),'_exp_fig.eps'), '-cmyk');
            export_fig(newFigure, strcat(handles.pathname,'plots', filesep,num2str(handles.eventselected),'_exp_fig.png'), '-cmyk');
        end
%         close(newFig)
    end
    
    % get the DNA level from the text box
    DNAdelta = str2num(get(handles.edit_dna_level,'String'));
    eventtrace_baseline = handles.rawtrace_vect(eventpoints(1):eventpoints(2)) - handles.event_info_temp(handles.eventselected,5);
    fold_count = 0;
    
    % normalize it so that each DNA level is an integer 
    event_levels_norm_raw = eventtrace_baseline/(-1.*DNAdelta);
    event_levels_norm_floor = floor(event_levels_norm_raw);
    event_levels_norm_rest = event_levels_norm_raw - event_levels_norm_floor;
    
    event_levels_norm = round(event_levels_norm_raw);
    
%     if event_levels_norm_rest(1) > 0.5
%         event_levels_norm(1) == event_levels_norm_floor(1) + 1;
%     else
%         event_levels_norm(1) == event_levels_norm_floor(1);
%     end
    
    event_levels_norm_full = zeros(length(event_levels_norm)+1,1);
    event_levels_norm_full(2:end,1) = event_levels_norm;
    num_levels = length(event_levels_norm);
    
    % get rounding limits
    round_up_limit = str2num(get(handles.edit_dna_round_up,'String')); % 
    round_down_limit = str2num(get(handles.edit_dna_round_down,'String')); % 
    
    % calculate the fold count, this part is obsolete, export to
    % OpenNanopore for level fitting
    for j=1:num_levels
        if event_levels_norm_rest(j) > round_up_limit % this value is sensitive to noise
            event_levels_norm_full(j+1) = event_levels_norm_floor(j) + 1;
        elseif event_levels_norm_rest(1) < round_down_limit % second value
            event_levels_norm_full(j+1) = event_levels_norm_floor(j);
        else
            event_levels_norm_full(j+1) = event_levels_norm_full(j);
        end
        % change in level
        current_diff = event_levels_norm_full(j+1)-event_levels_norm_full(j);
        if current_diff >= 1
            % count the change if larger than 1
            fold_count = fold_count + current_diff; 
        end
        
    end
    fold_count = fold_count-1;
    set(handles.edit_num_folds,'String',num2str(fold_count));
    
    timetemp_levels = timetemp;
    end_number = length(timetemp_levels);
    timetemp_levels(end_number+1)= timetemp_levels(end_number)+handles.timestep*1000;
    
    % this will plot the calculated levels
    if get(handles.check_show_levels,'Value')
        if get(handles.check_remove_baseline,'Value')
            if get(handles.check_use_extra_points,'Value') % extra points
                if get(handles.check_conductance,'Value') % conductance
                    plot(timetemp(handles.extrapoints+1:end-handles.extrapoints+1),((-1.*DNAdelta)*event_levels_norm_full(1:end))./voltage,'-r');
                else % no conductance
                    plot(timetemp(handles.extrapoints+1:end-handles.extrapoints+1),(-1.*DNAdelta)*event_levels_norm_full(1:end),'-r');
                end

            else % no extra points
                if get(handles.check_conductance,'Value') % conductance
                    plot(timetemp_levels,(-1.*DNAdelta)*event_levels_norm_full(1:end)./voltage,'-r');
                else % no conductance
                    plot(timetemp_levels,(-1.*DNAdelta)*event_levels_norm_full(1:end),'-r');
                end
            end
        else % baseline not removed
            if get(handles.check_use_extra_points,'Value') % extra points
                if get(handles.check_conductance,'Value') % conductance (baseline is removed)
                    plot(timetemp(handles.extrapoints+1:end-handles.extrapoints+1),((-1.*DNAdelta)*event_levels_norm_full(1:end))./voltage,'-r');
                else % no conductance
                    plot(timetemp(handles.extrapoints+1:end-handles.extrapoints+1),(-1.*DNAdelta)*event_levels_norm_full(1:end)+ handles.event_info_temp(handles.eventselected,5),'-r');
                end

            else % no extra points
                if get(handles.check_conductance,'Value') % conductance (baseline is removed)
                    plot(timetemp_levels,(-1.*DNAdelta)*event_levels_norm_full(1:end)./voltage,'-r');
                else % no conductance
                    plot(timetemp_levels,(-1.*DNAdelta)*event_levels_norm_full(1:end) + handles.event_info_temp(handles.eventselected,5),'-r');
                end
            end
        end
    end
    
    guidata(gcbo, handles); % save handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function listbox_eventsbox_CreateFcn(hObject, eventdata, handles)
    % Hint: listbox controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save currently selected events in standard format
function button_save_events_DB_Callback(hObject, ~, handles)

    % get output directory
    handles.outpathname = uigetdir; 

    set(handles.text_status,'String','Preparing matrices for output...parsing rawtrace...')

    % TODO this needs to be checked and updated
    indpos = 1;
    handles.extrapoints = str2num(get(handles.edit_extra_points,'String'));

    % pre allocate pointers for speed
    newstart = zeros(size(handles.event_info_temp,1),1);
    newend = zeros(size(handles.event_info_temp,1),1);

    % save sorted rawtrace
    for i=1:size(handles.event_info_temp,1) 
        % event start stop
        eventpoints = [(handles.event_info_temp(i,8)) (handles.event_info_temp(i,9))]; 
        % event duration in points
        eventlength = handles.event_info_temp(i,9) - handles.event_info_temp(i,8) + handles.extrapoints*2 + 1; 
        if isfield(handles,'rawtrace_vect') % only do this if rawtrace exists
            % create trace
            eventtrace = handles.rawtrace_vect(eventpoints(1)-handles.extrapoints:eventpoints(2)+handles.extrapoints); 
            % concatenate individual events together
            newtrace(indpos:indpos+eventlength-1) = eventtrace; 
            % vector with baseline removed
            newtrace_baseline(indpos:indpos+eventlength-1) = eventtrace - handles.event_info_temp(i,5); 
            if isfield(handles,'rawunfilteredtrace_vect') % do it for unfiltered data too?
                % create trace
                eventunfilteredtrace = handles.rawunfilteredtrace_vect(eventpoints(1)-handles.extrapoints:eventpoints(2)+handles.extrapoints); 
                % concatenate individual events together
                newunfilteredtrace(indpos:indpos+eventlength-1) = eventunfilteredtrace; 
            end
        else
            set(handles.text_status,'String','Rawtrace file not loaded!!')
        end
        % new start position
        newstart(i) = indpos + handles.extrapoints; 
        % new end position
        newend(i) = newstart(i) + handles.event_info_temp(i,9) - handles.event_info_temp(i,8); 
        % current position in trace in points
        indpos = indpos + eventlength; 
    end
    set(handles.text_status,'String','Preparing matrices for output...parsing rawtrace...done...writing...')

    % duplicate event info vector before changing it
    temp_out = handles.event_info_temp; 
    % update event start points
    temp_out(:,8) = newstart; 
    % update event end points
    temp_out(:,9) = newend; 

    % write events file
    dlmwrite(strcat(handles.outpathname, filesep,'analysis_events'), temp_out, 'delimiter', '\t', ...
             'precision', '%.8f'); 
    if isfield(handles,'events_minmax_temp') % write minmax file
        dlmwrite(strcat(handles.outpathname, filesep,'analysis_events_min_max'), handles.events_minmax_temp, 'delimiter', '\t', ...
                 'precision', '%.8f'); 
    end
    % write extra info file
    dlmwrite(strcat(handles.outpathname, filesep,'analysis_extra_info.txt'), handles.event_more_info_temp, 'delimiter', '\t', ...
             'precision', '%.8f'); 
    if isfield(handles,'rawtrace_vect')
        % write rawtrace file
        dlmwrite(strcat(handles.outpathname, filesep,'analysis_rawtrace'), newtrace, 'delimiter', '\n', ...
             'precision', '%.8f'); 
        if isfield(handles,'rawunfilteredtrace_vect')
            % write unfiltered rawtrace
            dlmwrite(strcat(handles.outpathname, filesep,'analysis_unfiltered_rawtrace'), newunfilteredtrace, 'delimiter', '\n', ...
             'precision', '%.8f'); 
        end
        % write OpenNanopore file
        newtrace_baseline = transpose(newtrace_baseline);
        save(strcat(handles.outpathname, filesep,'OpenNanoporeEvents.mat'), 'newtrace_baseline');
    end

    % copy ~detect.par
    copyfile(strcat(handles.pathname, filesep,'~detect.par'),strcat(handles.outpathname, filesep,'~detect.par'))
    % copy legacy file
    copyfile(strcat(handles.pathname, filesep,'analysis_parameters.txt'),strcat(handles.outpathname, filesep,'analysis_parameters.txt'))

    set(handles.text_status,'String','Preparing matrices for output...parsing rawtrace...done...writing...done')
    guidata(hObject, handles); % update handles



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in an events database
function button_load_events_file_Callback(hObject, ~, handles)

    clear rawtrace_vect

    set(handles.text_status,'String','Waiting for file...')
    % get the events file
    [handles.filename, handles.pathname] = uigetfile( ...
    {  '*',  'All Files (*.*)'; ...
       '*_events','Events File'}, ...
       'Pick a file'); 
    % legacy code, before standardized naming
    [fileactual filediscard] = strread(handles.filename, '%s %s', 'delimiter','_'); 

    % determine rawtrace file name
    handles.filenametrace = char(strcat(fileactual,'_rawtrace')); 
    % in case old selection value higher than new total (which causes events box to disappear)
    set(handles.listbox_eventsbox,'Value',1); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.text_status,'String','Opening events file...')
    % open events file
    fid = fopen(strcat(handles.pathname,handles.filename)); 
    if fid == -1
        disp(strcat(handles.filename,' could not be loaded!!')) % not good
    else
        % grab the event info
        handles.event_info = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f')); 
        % contents of this vector (for each event):
    %                 1 amplitude_minima_maxima
    %                 2 maximum	Amplitude
    %                 3 integral
    %                 4 FWHM_dwell
    %                 5 baseline
    %                 6 detectionlevel
    %                 7 eventtype UpDown	down =0 up=1
    %                 8 startpoint
    %                 9 endpoint
    end
    fclose(fid);
    set(handles.text_status,'String','Opening events file...done...opening rawtrace file...')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read analysis_extra_info.txt
    % 
    filenameextrainfo = strcat(handles.pathname, filesep,'analysis_extra_info.txt');
    if isequal(exist(filenameextrainfo),2) % if the extra info file exists
        % open it
        fid = fopen(filenameextrainfo); 
        if fid == -1
            disp(strcat(filenameextrainfo,' could not be loaded!!'))
            % can't do event rate calculations, turn GUI_eventrate off
            set(handles.button_event_rate_mgmt,'Enable','off') 
        else
            set(handles.button_event_rate_mgmt,'Enable','on')
            handles.event_more_info = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f')); % get the extra info
            handles.event_more_info_temp = handles.event_more_info; % select everything
            % this contains :
            % 1 file number
            % 2 resistance for this segment
            % 3 time - local within segment
            % 4 time - global from start of experiment (first segment in directory)
            % 5 number of local min
            % 6 index first local min
            % 7 index last local min
            % 8 dwellold (between start and stop)
            % 9 amplitude good
        end

        fclose(fid);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this is commented out because it's not used anywhere
%     filenameeventsource = strcat(handles.pathname, filesep,'analysis_files_info.txt');
%     if isequal(exist(filenameeventsource),2) % if it exists
%         [handles.trace_file_names_array,~]= readtext(filenameeventsource, '\t', '', '"', '');
%         % handles.trace_file_names_array is an array of the name of the
%         % file where every event is from
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read in analysis_events_min_max
    filenameminmaxinfo = strcat(handles.pathname, filesep,'analysis_events_min_max');
    if isequal(exist(filenameminmaxinfo),2) % if the extra info file exists
        % open it
        fid = fopen(filenameminmaxinfo); 
        if fid == -1 % catch problems
            disp(strcat(filenameminmaxinfo,' could not be loaded!!'))
        else
            % get the data
            handles.events_minmax = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f %f')); 
            % this is for backwards compatibility
            handles.events_minmax(:,8) = abs(handles.events_minmax(:,8)); 
            % this is for backwards compatibility
            handles.events_minmax(:,10) = abs(handles.events_minmax(:,10)); 
            % select everything
            handles.events_minmax_temp = handles.events_minmax; 
        end
        %this vector has:
        %1 min inside event
        %2 max inside event
        %3 min in extrapoints
        %4 max in extrapoints
        %5 std in extrapoints
        %6 mean in extrapoints
        %7 std in extrapoints start
        %8 mean extrapoints start
        %9 std in extrapoints end
        %10 mean extrapoints end

        fclose(fid);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read in analysis_parameters.txt
    filenameparam = strcat(handles.pathname, filesep,'analysis_parameters.txt');
    if isequal(exist(filenameparam),2) % if the paramter file exists
        fid = fopen(filenameparam);
        if fid == -1 % catch problem
            disp(strcat(filenameparam,' could not be loaded!!'))
        else % file loaded successfully
            line = fgetl(fid); % skip Transalyzer Analysis
            line = fgetl(fid); % skip version
            line = fgetl(fid); % skip Analysis Date
            line = fgetl(fid); % date, not used anywhere
            line = fgetl(fid); % skip
            handles.voltage = str2double(fgetl(fid)); % Voltage (mV)
            set(handles.edit_voltage,'String',num2str(handles.voltage));
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
            set(handles.edit_extra_points,'String',num2str(handles.readextrapoints));
            line = fgetl(fid); % skip
            handles.totaleventsfound = str2double(fgetl(fid)); % Total Events Found
            line = fgetl(fid); % skip
            handles.totaltimetrace = str2double(fgetl(fid)); % Total Time
        end
        fclose(fid); % close the file
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read in unfiltered current trace file (rawtrace)
    if get(handles.check_load_rawtrace,'Value') % only load rawtrace if user wants to
        filenameunfilttrace = strcat(handles.pathname, filesep,'analysis_unfiltered_rawtrace');
        if isequal(exist(filenameunfilttrace),2) % check file exists
            fid = fopen(filenameunfilttrace); % open it
            if fid == -1 % catch problems
                disp(strcat(filenameunfilttrace,' could not be loaded!!'))
                % need unfiltered data for this to work
                set(handles.check_overlay_unfiltered,'Enable','off') 
                % need unfiltered data for this to work
                set(handles.button_ghost,'Enable','off')
            else
                set(handles.check_overlay_unfiltered,'Enable','on')
                set(handles.button_ghost,'Enable','on')
                handles.rawunfilteredtrace_vect = cell2mat(textscan(fid, '%f')); % read it in
            end
            fclose(fid);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % The old comma to period conversion code was here.
    % Removed & all computers switched to using period
    % for the decimal separator.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if get(handles.check_multiply_by_ten,'Value') == 1
        % multiply by 10, this is legacy support for the old analysis code
        handles.event_info = [handles.event_info(:,1).*10 handles.event_info(:,2).*10 handles.event_info(:,3).*10 handles.event_info(:,4) handles.event_info(:,5).*10 handles.event_info(:,6).*10 handles.event_info(:,7) handles.event_info(:,8) handles.event_info(:,9)];
    end

    % handles.event_info contains everything, this is always used at start of
    % sorting
    % handles.event_info_temp contains sorted data
    
    % select everything to begin
    handles.event_info_temp = handles.event_info; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read in current trace file (rawtrace)
    if get(handles.check_load_rawtrace,'Value') % if user want to load rawtrace
        fid = fopen(strcat(handles.pathname,handles.filenametrace));
        if fid == -1 % catch problems
                disp(strcat(handles.filename,' could not be loaded!!'))
        else
            % read in rawtrace file
            handles.rawtrace_vect = cell2mat(textscan(fid, '%f')); 
        end
        if get(handles.check_multiply_by_ten,'Value') == 1
            % legacy support
            handles.rawtrace_vect = handles.rawtrace_vect.*10; 
        end
        fclose(fid);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update vatious text boxes / radio buttons with the new data, and/or default values and
    % states
    set(handles.text_status,'String','Opening events file...done...opening rawtrace file...done...plotting')
    set(handles.text_file_box,'String',strcat(handles.pathname,handles.filename))
    set(handles.listbox_eventsbox,'String',num2str(handles.event_info_temp))
    set(handles.text_number_of_events,'String',num2str(size(handles.event_info_temp,1)))
    set(handles.radio_2_scatter,'Value',1);
    set(handles.radio_2_base,'Value',0);
    set(handles.radio_2_dwell,'Value',0);
    set(handles.radio_2_amplitude,'Value',0);
    set(handles.radio_2_current,'Value',0);
    set(handles.radio_2_scatter_integral,'Value',0);
    set(handles.radio_2_integral_hist,'Value',0);

    % plot scatter after loading
    fun_plot_scatter(handles,'X') 

    % set min and max Dwell values in the event selection criteria boxes
    switch get(handles.pop_time_type,'Value') 
        case 1 % FWHM_dwell
            set(handles.edit_select_dwell_min,'String',num2str(min(handles.event_info_temp(:,4))))
            set(handles.edit_select_dwell_max,'String',num2str(max(handles.event_info_temp(:,4))))
        case 2 % dwellold (between start and stop)
            set(handles.edit_select_dwell_min,'String',num2str(min(handles.event_more_info_temp(:,8))))
            set(handles.edit_select_dwell_max,'String',num2str(max(handles.event_more_info_temp(:,8))))
        otherwise
    end

    % set min and max Amplitude values in the event selection criteria boxes
    switch get(handles.pop_amplitude_type,'Value') 
        case 1 % amplitude_minima_maxima
            set(handles.edit_select_amplitude_min,'String',num2str(min(handles.event_info_temp(:,1))))
            set(handles.edit_select_amplitude_max,'String',num2str(max(handles.event_info_temp(:,1))))
        case 2 % Max Amplitude
            set(handles.edit_select_amplitude_min,'String',num2str(min(handles.event_info_temp(:,2))))
            set(handles.edit_select_amplitude_max,'String',num2str(max(handles.event_info_temp(:,2))))
        case 3 % Integral / FWHM
            set(handles.edit_select_amplitude_min,'String',num2str(min(handles.event_more_info_temp(:,9))))
            set(handles.edit_select_amplitude_max,'String',num2str(max(handles.event_more_info_temp(:,9))))
        otherwise
    end

    % set min and max Integral values in the event selection criteria boxes
    set(handles.edit_select_integral_min,'String',num2str(min(handles.event_info_temp(:,3))))
    set(handles.edit_select_integral_max,'String',num2str(max(handles.event_info_temp(:,3))))

    % set min and max Baseline values in the event selection criteria boxes
    set(handles.edit_select_baseline_min,'String',num2str(min(handles.event_info_temp(:,5))))
    set(handles.edit_select_baseline_max,'String',num2str(max(handles.event_info_temp(:,5))))

    % set min and max Event # in the event selection criteria boxes
    set(handles.edit_select_event_number_min,'String','1')
    set(handles.edit_select_event_number_max,'String',num2str(length(handles.event_info(:,5)))) % last event #

    if isfield(handles,'events_minmax')
        % set min and max current points in the event selection criteria boxes
        set(handles.edit_select_firstextrapointmeanmax,'String',num2str(max(handles.events_minmax(:,8))))
        set(handles.edit_select_lastextrapointmeanmax,'String',num2str(max(handles.events_minmax(:,10))))
    end

    set(handles.text_status,'String','Opening events file...done...opening rawtrace file...done...plotting...done')
    fun_update_detection_level(handles);
    guidata(hObject, handles); % update handles



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_exit_Callback(hObject, eventdata, handles)

    delete(gcf) % close it up


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_scatter(handles,location)
    %shading faceted
    %grid on
    
    if strcmp(location,'X') % plot in GUI?
        cla(handles.axes_graphs,'reset')
        axes(handles.axes_graphs)
    end
    if strcmp(location,'n') % make new figure?
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        figure()
    end
    % read in voltage
    voltage = str2double(get(handles.edit_voltage,'String'))/1000; 
    
    amplitude = fun_get_amplitude(handles);

    switch get(handles.pop_time_type,'Value') % type of Dwell calc to use?
        case 1 % FWHM_dwell
            dwelldata = handles.event_info_temp(:,4);
        case 2 % dwellold (between start and stop)
            dwelldata = handles.event_more_info_temp(:,8);
        otherwise
    end

    if strcmp(location,'h') % create a scatter hist plot, called by mass export button press
        bins = str2num(get(handles.edit_2_numbins,'String'));
        pp3 = scatterhist(dwelldata,amplitude,[bins bins]);
    else % otherwise use normal scatter plot
        pp3 = scatter(dwelldata,amplitude,'or');
    end
    if get(handles.check_dwell_X_log,'Value') % should we use log scale?
        set( gca                       , ...
        'XScale'   , 'log' );
    end
    
    
    % add axis labels
    hXLabel = xlabel('Time (ms)');
    if get(handles.check_relative,'Value')
        hYLabel = ylabel('Relative Blockade');
    elseif get(handles.check_conductance,'Value')
        hYLabel = ylabel('Conductance (nS)');
    else
        hYLabel = ylabel('Current Blockade (nA)');
    end

    if strcmp(location,'n') % use limits same as on old plot in GUI
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
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
      'XGrid'       , 'on'      , ...
      'XColor'      , [0 0 0], ...
      'YColor'      , [0 0 0], ...
      'LineWidth'   , 1         );

    if strcmp(location,'s') % save the data
        dir_str = strcat(handles.pathname,'plots', filesep); % in the plots dir
        if(~exist(dir_str, 'dir'))
            mkdir(dir_str);
        end
        file_str = strcat(dir_str,'plot_data_scatter.dat'); % save dat file with scatter data
        dataoutput = [dwelldata, amplitude];
        dlmwrite(file_str, dataoutput, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    end

    if strcmp(location,'X') 
        % if mouse click on plot in GUI, get mouse data
        set(pp3, 'ButtonDownFcn', {@plot2_ButtonDownFcn, handles}); 

    end
    
    % set your plot limits
    ylimits = get(handles.axes_graphs,'YLim');
    xlimits = get(handles.axes_graphs,'XLim');
    
    % update text boxes for limits
    set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
    set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
    set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
    set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
    
    % TODO ugly, using globals to pass data
    global scatterXlim; 
    scatterXlim = xlimits;
    global scatterYlim; 
    scatterYlim= ylimits;


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun_plot_scatter_heat(handles,location)
    global scatterXlim; % TODO ugly, using globals
    global scatterYlim; % TODO ugly, using globals
    
    if strcmp(location,'X') % plot in GUI?
        cla(handles.axes_graphs,'reset')
        axes(handles.axes_graphs)
    end
    
    ylimits = scatterYlim;
    xlimits = scatterXlim;
    
    if strcmp(location,'n') % make new figure?
        figure()
    end
    
    % read in voltage
    voltage = str2double(get(handles.edit_voltage,'String'))/1000; 

    amplitude = fun_get_amplitude(handles);

    switch get(handles.pop_time_type,'Value') % type of Dwell calc to use?
        case 1 % FWHM_dwell
            dwelldata = handles.event_info_temp(:,4);
        case 2 % dwellold (between start and stop)
            dwelldata = handles.event_more_info_temp(:,8);
        otherwise
    end

    % number of bins
    bins = str2num(get(handles.edit_heatbins,'String'));
    
    if strcmp(location,'h') % scatter hist plot, called by mass export button press
        
        pp3 = scatterhist(dwelldata,amplitude,[bins bins]);
        
    else % otherwise use normal scatter heat plot

        xi = linspace(min(dwelldata(:)),max(dwelldata(:)),bins);
        yi = linspace(min(amplitude(:)),max(amplitude(:)),bins);

        xr = interp1(xi,1:numel(xi),dwelldata,'nearest')';
        yr = interp1(yi,1:numel(yi),amplitude,'nearest')';

        scatter2dheat = log10(accumarray([xr' yr'], 1, [bins bins]));
        scatter2dheat(isinf(scatter2dheat))=0; % replace inf with 0

        pp3 = surf(xi,yi,scatter2dheat');
        
        colorbar
        %this is a nice colormap
        colormap([1 1 1;0.949999988079071 1 1;0.899999976158142 1 1;0.850000023841858 1 1;0.800000011920929 1 1;0.75 1 1;0.699999988079071 1 1;0.649999976158142 1 1;0.600000023841858 1 1;0.550000011920929 1 1;0.5 1 1;0.449999988079071 1 1;0.400000005960464 1 1;0.349999994039536 1 1;0.300000011920929 1 1;0.25 1 1;0.200000002980232 1 1;0.150000005960464 1 1;0.100000001490116 1 1;0.0500000007450581 1 1;0 1 1;0.0476190485060215 1 0.952380955219269;0.095238097012043 1 0.904761910438538;0.142857149243355 1 0.857142865657806;0.190476194024086 1 0.809523820877075;0.238095238804817 1 0.761904776096344;0.28571429848671 1 0.714285731315613;0.333333343267441 1 0.666666686534882;0.380952388048172 1 0.61904764175415;0.428571432828903 1 0.571428596973419;0.476190477609634 1 0.523809552192688;0.523809552192688 1 0.476190477609634;0.571428596973419 1 0.428571432828903;0.61904764175415 1 0.380952388048172;0.666666686534882 1 0.333333343267441;0.714285731315613 1 0.28571429848671;0.761904776096344 1 0.238095238804817;0.809523820877075 1 0.190476194024086;0.857142865657806 1 0.142857149243355;0.904761910438538 1 0.095238097012043;0.952380955219269 1 0.0476190485060215;1 1 0;1 0.954545438289642 0;1 0.909090936183929 0;1 0.863636374473572 0;1 0.818181812763214 0;1 0.772727251052856 0;1 0.727272748947144 0;1 0.681818187236786 0;1 0.636363625526428 0;1 0.590909063816071 0;1 0.545454561710358 0;1 0.5 0;1 0.454545468091965 0;1 0.409090906381607 0;1 0.363636374473572 0;1 0.318181812763214 0;1 0.272727280855179 0;1 0.227272734045982 0;1 0.181818187236786 0;1 0.136363640427589 0;1 0.0909090936183929 0;1 0.0454545468091965 0;1 0 0]);
 
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2) 0 max(max(scatter2dheat))]);
        view(2)
        grid off

        shading flat
        %shading interp

    end
    if get(handles.check_dwell_X_log,'Value') % use log axis?
        set( gca                       , ...
        'XScale'   , 'log' );
    end

    % add axis labels
    hXLabel = xlabel('Time (ms)');
    if get(handles.check_conductance,'Value')
        hYLabel = ylabel('Conductance (nS)');
    elseif get(handles.check_relative,'Value')
        hYLabel = ylabel('Relative Blockade');
    else
        hYLabel = ylabel('Current Blockade (nA)');
    end

    if strcmp(location,'n') % use limits same as on old plot in GUI
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
    end

    % style
    set( gca                       , ...
        'FontName'   , 'Arial' );
    set([hXLabel, hYLabel], ...
        'FontName'   , 'Arial');
    set([hXLabel, hYLabel]  , ...
        'FontSize'   , 10          );

    set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'on'      , ...
      'YMinorTick'  , 'on'      , ...
      'XColor'      , [0 0 0], ...
      'YColor'      , [0 0 0], ...
      'LineWidth'   , 1         );

    if strcmp(location,'s') % save plot data?
        dir_str = strcat(handles.pathname,'plots', filesep); % in the plots dir
        if(~exist(dir_str, 'dir'))
            mkdir(dir_str);
        end
        file_str = strcat(dir_str,'plot_data_scatter.dat'); % save dat file with scatter data
        dataoutput = [dwelldata, amplitude];
        dlmwrite(file_str, dataoutput, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    end
    if strcmp(location,'n') % save figure?
        dir_str = strcat(handles.pathname,'plots', filesep); % in the plots dir
        if(~exist(dir_str, 'dir'))
            mkdir(dir_str);
        end
        if exist('export_fig') == 2
            set(gcf, 'Color', 'w');
            export_fig(gcf, strcat(handles.pathname,'plots', filesep,'plot_scatter_heat_exp_fig.eps'), '-cmyk');
            export_fig(gcf, strcat(handles.pathname,'plots', filesep,'plot_scatter_heat_exp_fig.png'), '-cmyk');
        end
    end
    if strcmp(location,'X') 
        % if mouse click on plot, save details
        set(pp3, 'ButtonDownFcn', {@plot2_ButtonDownFcn, handles}); 

        % set your plot limits
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
        set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
        set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
        set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
    end 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_dwell_min_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_dwell_min_CreateFcn(hObject, eventdata, handles)

    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_dwell_max_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_dwell_max_CreateFcn(hObject, eventdata, handles)

    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_amplitude_min_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_amplitude_min_CreateFcn(hObject, eventdata, handles)

    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_amplitude_max_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_amplitude_max_CreateFcn(hObject, eventdata, handles)

    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_baseline_min_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_baseline_min_CreateFcn(hObject, eventdata, handles)

    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_baseline_max_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_baseline_max_CreateFcn(hObject, eventdata, handles)

    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this selects events based on the limits entered into the text boxes
function button_select_events_Callback(hObject, ~, handles)

    set(handles.text_status,'String','Sorting Events...wait')

    selectbasemin = str2num(get(handles.edit_select_baseline_min,'String')); % min baseline
    selectbasemax = str2num(get(handles.edit_select_baseline_max,'String')); % max baseline
    selectevenummin = str2num(get(handles.edit_select_event_number_min,'String')); % min event #
    selectevenummax = str2num(get(handles.edit_select_event_number_max,'String')); % max event #
    selectintmin = str2num(get(handles.edit_select_integral_min,'String')); % min integral
    selectintmax = str2num(get(handles.edit_select_integral_max,'String')); % max intrgeal
    selectfirstexptmenmax = str2num(get(handles.edit_select_firstextrapointmeanmax,'String')); % max first extrapoint mean
    selectlastexptmenmax = str2num(get(handles.edit_select_lastextrapointmeanmax,'String')); % max last extrapoint mean

    % selection area is a box in scatter plot
    if get(handles.radio_box_mode,'Value') == 1 

        select_dwell_min = str2num(get(handles.edit_select_dwell_min,'String')); % min dwell
        select_dwell_max = str2num(get(handles.edit_select_dwell_max,'String')); % max dwell

        select_amp_min = str2num(get(handles.edit_select_amplitude_min,'String')); % min amp
        select_amp_max = str2num(get(handles.edit_select_amplitude_max,'String')); % max amp

        % 6 cases in total 3x amplitude X 2x dwell time
        % in each case apply selection on proper variables

        switch get(handles.pop_amplitude_type,'Value') 
            case 1 % amplitude_minima_maxima
                switch get(handles.pop_time_type,'Value') 
                    case 1 % FWHM_dwell
                        selectindex = find(handles.event_info(:,4) >= select_dwell_min & handles.event_info(:,4) <= select_dwell_max &...
                            handles.event_info(:,1) >= select_amp_min & handles.event_info(:,1) <= select_amp_max &...
                            handles.event_info(:,5) >= selectbasemin & handles.event_info(:,5) <= selectbasemax &...
                            handles.event_info(:,3) >= selectintmin & handles.event_info(:,3) <= selectintmax &...
                            handles.events_minmax(:,8) <= selectfirstexptmenmax & handles.events_minmax(:,10) <= selectlastexptmenmax);
                    case 2 % dwellold (between start and stop)
                        selectindex = find(handles.event_more_info(:,8) >= select_dwell_min & handles.event_more_info(:,8) <= select_dwell_max &...
                            handles.event_info(:,1) >= select_amp_min & handles.event_info(:,1) <= select_amp_max &...
                            handles.event_info(:,5) >= selectbasemin & handles.event_info(:,5) <= selectbasemax &...
                            handles.event_info(:,3) >= selectintmin & handles.event_info(:,3) <= selectintmax &...
                            handles.events_minmax(:,8) <= selectfirstexptmenmax & handles.events_minmax(:,10) <= selectlastexptmenmax);
                    otherwise
                end
            case 2 % Max Amplitude
                switch get(handles.pop_time_type,'Value') 
                    case 1 % FWHM_dwell
                        selectindex = find(handles.event_info(:,4) >= select_dwell_min & handles.event_info(:,4) <= select_dwell_max &...
                            handles.event_info(:,2) >= select_amp_min & handles.event_info(:,2) <= select_amp_max &...
                            handles.event_info(:,5) >= selectbasemin & handles.event_info(:,5) <= selectbasemax &...
                            handles.event_info(:,3) >= selectintmin & handles.event_info(:,3) <= selectintmax &...
                            handles.events_minmax(:,8) <= selectfirstexptmenmax & handles.events_minmax(:,10) <= selectlastexptmenmax);         
                    case 2 % dwellold (between start and stop)
                        selectindex = find(handles.event_more_info(:,8) >= select_dwell_min & handles.event_more_info(:,8) <= select_dwell_max &...
                            handles.event_info(:,2) >= select_amp_min & handles.event_info(:,2) <= select_amp_max &...
                            handles.event_info(:,5) >= selectbasemin & handles.event_info(:,5) <= selectbasemax &...
                            handles.event_info(:,3) >= selectintmin & handles.event_info(:,3) <= selectintmax &...
                            handles.events_minmax(:,8) <= selectfirstexptmenmax & handles.events_minmax(:,10) <= selectlastexptmenmax);         
                    otherwise
                end
            case 3 % Integral / FWHM
                switch get(handles.pop_time_type,'Value') 
                    case 1 % FWHM_dwell % this is the default selection
                        selectindex = find(handles.event_info(:,4) >= select_dwell_min & handles.event_info(:,4) <= select_dwell_max &...
                            handles.event_more_info(:,9) >= select_amp_min & handles.event_more_info(:,9) <= select_amp_max &...
                            handles.event_info(:,5) >= selectbasemin & handles.event_info(:,5) <= selectbasemax &...
                            handles.event_info(:,3) >= selectintmin & handles.event_info(:,3) <= selectintmax &...
                            handles.events_minmax(:,8) <= selectfirstexptmenmax & handles.events_minmax(:,10) <= selectlastexptmenmax);         
                    case 2 % dwellold (between start and stop)
                        selectindex = find(handles.event_more_info(:,8) >= select_dwell_min & handles.event_more_info(:,8) <= select_dwell_max &...
                            handles.event_more_info(:,9) >= select_amp_min & handles.event_more_info(:,9) <= select_amp_max &...
                            handles.event_info(:,5) >= selectbasemin & handles.event_info(:,5) <= selectbasemax &...
                            handles.event_info(:,3) >= selectintmin & handles.event_info(:,3) <= selectintmax &...
                            handles.events_minmax(:,8) <= selectfirstexptmenmax & handles.events_minmax(:,10) <= selectlastexptmenmax);         
                    otherwise
                end
            otherwise
        end
    end

    % selection area is ellipse in scatter plot
    % TODO this seems to have some issues, check
    if get(handles.radio_ellipse_mode,'Value') == 1 

        % remember for an ellipse:
        % ((x-xc)*cos(t)-(y-yc)*sin(t)).^2/a^2 + ...
        % ((x-xc)*sin(t)+(y-yc)*cos(t)).^2/b^2 <= 1

        mindwell = str2num(get(handles.edit_select_dwell_min,'String')); % dwell select min
        maxdwell = str2num(get(handles.edit_select_dwell_max,'String')); % dwell select max

        minavg = str2num(get(handles.edit_select_amplitude_min,'String')); % avg select min
        maxavg = str2num(get(handles.edit_select_amplitude_max,'String')); % avg select max

        widthdwell = maxdwell - mindwell; % ellipse width
        heightavg = maxavg - minavg; % ellipse height

        % 6 cases in total 3x amplitude X 2x dwell time
        % in each case apply selection on proper variables

        switch get(handles.pop_amplitude_type,'Value') 
            case 1 % amplitude_minima_maxima (between local minima)
                switch get(handles.pop_time_type,'Value') 
                    case 1 % FWHM_dwell
                        selectindex = find((handles.event_info(:,4)-(widthdwell/2+mindwell)).^2./(widthdwell/2)^2 +...
                            (handles.event_info(:,1)-(heightavg/2+minavg)).^2./(heightavg/2)^2 <=  1 &...
                            handles.event_info(:,5) > selectbasemin & handles.event_info(:,5) < selectbasemax &...
                            handles.event_info(:,3) >= selectintmin & handles.event_info(:,3) <= selectintmax &...
                            handles.events_minmax(:,8) <= selectfirstexptmenmax & handles.events_minmax(:,10) <= selectlastexptmenmax);
                    case 2 % dwellold (between start and stop)
                        selectindex = find((handles.event_more_info(:,8)-(widthdwell/2+mindwell)).^2./(widthdwell/2)^2 +...
                            (handles.event_info(:,1)-(heightavg/2+minavg)).^2./(heightavg/2)^2 <=  1 &...
                            handles.event_info(:,5) > selectbasemin & handles.event_info(:,5) < selectbasemax &...
                            handles.event_info(:,3) >= selectintmin & handles.event_info(:,3) <= selectintmax &...
                            handles.events_minmax(:,8) <= selectfirstexptmenmax & handles.events_minmax(:,10) <= selectlastexptmenmax);
                    otherwise
                end
            case 2 % Max Amplitude (for proteins and short DNA)
                switch get(handles.pop_time_type,'Value') 
                    case 1 % FWHM_dwell
                        selectindex = find((handles.event_info(:,4)-(widthdwell/2+mindwell)).^2./(widthdwell/2)^2 +...
                            (handles.event_info(:,2)-(heightavg/2+minavg)).^2./(heightavg/2)^2 <=  1 &...
                            handles.event_info(:,5) > selectbasemin & handles.event_info(:,5) < selectbasemax &...
                            handles.event_info(:,3) >= selectintmin & handles.event_info(:,3) <= selectintmax &...
                            handles.events_minmax(:,8) <= selectfirstexptmenmax & handles.events_minmax(:,10) <= selectlastexptmenmax);
                    case 2 % dwellold (between start and stop)
                        selectindex = find((handles.event_more_info(:,8)-(widthdwell/2+mindwell)).^2./(widthdwell/2)^2 +...
                            (handles.event_info(:,2)-(heightavg/2+minavg)).^2./(heightavg/2)^2 <=  1 &...
                            handles.event_info(:,5) > selectbasemin & handles.event_info(:,5) < selectbasemax &...
                            handles.event_info(:,3) >= selectintmin & handles.event_info(:,3) <= selectintmax &...
                            handles.events_minmax(:,8) <= selectfirstexptmenmax & handles.events_minmax(:,10) <= selectlastexptmenmax);
                    otherwise
                end

            case 3 % Integral / FWHM
                switch get(handles.pop_time_type,'Value') 
                    case 1 % FWHM_dwell
                        selectindex = find((handles.event_info(:,4)-(widthdwell/2+mindwell)).^2./(widthdwell/2)^2 +...
                            (handles.event_more_info(:,9)-(heightavg/2+minavg)).^2./(heightavg/2)^2 <=  1 &...
                            handles.event_info(:,5) > selectbasemin & handles.event_info(:,5) < selectbasemax &...
                            handles.event_info(:,3) >= selectintmin & handles.event_info(:,3) <= selectintmax &...
                            handles.events_minmax(:,8) <= selectfirstexptmenmax & handles.events_minmax(:,10) <= selectlastexptmenmax);
                    case 2 % dwellold (between start and stop)
                        selectindex = find((handles.event_more_info(:,8)-(widthdwell/2+mindwell)).^2./(widthdwell/2)^2 +...
                            (handles.event_more_info(:,9)-(heightavg/2+minavg)).^2./(heightavg/2)^2 <=  1 &...
                            handles.event_info(:,5) > selectbasemin & handles.event_info(:,5) < selectbasemax &...
                            handles.event_info(:,3) >= selectintmin & handles.event_info(:,3) <= selectintmax &...
                            handles.events_minmax(:,8) <= selectfirstexptmenmax & handles.events_minmax(:,10) <= selectlastexptmenmax);
                    otherwise
                end
            otherwise
        end
    end

    % event # selection is done here
    % find events between min and max event numbers and get position of the
    % events
    selectindex = selectindex(selectindex >= selectevenummin & selectindex <= selectevenummax);

    % event # selection is completed here
    handles.event_info_temp = handles.event_info(selectindex,:);
    if isfield(handles,'event_more_info')
        handles.event_more_info_temp = handles.event_more_info(selectindex,:); % select other parameters too
    end

    if isfield(handles,'events_minmax')
        handles.events_minmax_temp = handles.events_minmax(selectindex,:);
    end

    % TODO this has not been enabled yet, old obsolete code, redo
    if get(handles.check_use_fold_counts,'Value') % max min current points done here

        folding_count_temp = handles.folding_count(selectindex);

        minmaxindex = find((folding_count_temp >= str2double(get(handles.edit_select_fold_min,'String')) & folding_count_temp <= str2double(get(handles.edit_select_fold_max,'String'))));

        handles.event_info_temp = handles.event_info_temp(minmaxindex,:);
        if isfield(handles,'event_more_info')
            handles.event_more_info_temp = handles.event_more_info_temp(minmaxindex,:); % select other parameters too
        end
        if isfield(handles,'events_minmax')
            handles.events_minmax_temp = handles.events_minmax_temp(minmaxindex,:);
        end
    end

    % set back to 1, or events box may disappear
    set(handles.listbox_eventsbox,'Value',1) 
    % update events box
    set(handles.listbox_eventsbox,'String',num2str(handles.event_info_temp)) 
    % update number of events
    set(handles.text_number_of_events,'String',num2str(size(handles.event_info_temp,1))) 
    set(handles.text_status,'String','Sorting Events...wait...done...plotting...')

    % update radio button states
    set(handles.radio_2_scatter,'Value',1);
    set(handles.radio_2_base,'Value',0);
    set(handles.radio_2_dwell,'Value',0);
    set(handles.radio_2_R_pore,'Value',0);
    set(handles.radio_2_amplitude,'Value',0);
    set(handles.radio_2_current,'Value',0);
    set(handles.radio_2_scatter_integral,'Value',0);
    set(handles.radio_2_integral_hist,'Value',0);

    fun_plot_scatter(handles,'X') % plot new histogram in GUI
    set(handles.text_status,'String','Sorting Events...wait...done...plotting...done')
    fun_update_detection_level(handles);
    guidata(hObject, handles); % update the handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_extra_points_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_extra_points_CreateFcn(hObject, eventdata, handles)

    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_plot2xmin_Callback(hObject, eventdata, handles) 

    % set the min X axis value

    set(handles.text_status,'String','Setting Plot2 Xmin')
    xlimits = get(handles.axes_graphs,'XLim');
    minvalue = str2double(get(handles.edit_plot2xmin,'String'));
    set(handles.axes_graphs,'XLim',[minvalue xlimits(2)]);
    guidata(hObject, handles); % update handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_plot2xmin_CreateFcn(hObject, eventdata, handles)

    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_plot2xmax_Callback(hObject, ~, handles) 

    % set the max X axis value

    set(handles.text_status,'String','Setting Plot2 Xmax')
    xlimits = get(handles.axes_graphs,'XLim');
    maxvalue = str2double(get(handles.edit_plot2xmax,'String'));
    set(handles.axes_graphs,'XLim',[xlimits(1) maxvalue]);
    guidata(hObject, handles); % update handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_plot2xmax_CreateFcn(hObject, eventdata, handles)

    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uipanel_plot_2_type_SelectionChangeFcn(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_box_mode_Callback(hObject, eventdata, handles)

    %set(handles.radio_ellipse_mode,'Value',0); % unselect other radio button
    %set(handles.radio_box_mode,'Value',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_ellipse_mode_Callback(hObject, eventdata, handles)

    %set(handles.radio_box_mode,'Value',0);  % unselect other radio button
    %set(handles.radio_ellipse_mode,'Value',1);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this will draw the current dwell, amplitude, and integral limits in the
% current plot
function button_draw_limits_Callback(hObject, ~, handles)


    voltage = str2double(get(handles.edit_voltage,'String'))/1000; % conductance case

    % check what type of plot this is
    if get(handles.radio_2_scatter,'Value') == 1 % draw the selection limits
        set(handles.text_status,'String','Drawing Scatter Limits')
        xlimits = get(handles.axes_graphs,'XLim'); % dwell
        ylimits = get(handles.axes_graphs,'YLim'); % amp
        axes(handles.axes_graphs) % axes_graphs is the destination

        fun_plot_scatter(handles,'X') % update scatter plot

        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
        set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
        set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
        set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
        set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))

        mindwell = str2num(get(handles.edit_select_dwell_min,'String')); % dwell select min
        maxdwell = str2num(get(handles.edit_select_dwell_max,'String')); % dwell select max
        minavg = str2num(get(handles.edit_select_amplitude_min,'String')); % amp select min
        maxavg = str2num(get(handles.edit_select_amplitude_max,'String')); % amp select max
        
        %%%%%%%%%%%%%%%%%%%%%%
        % draw integral limits
        minint = str2num(get(handles.edit_select_integral_min,'String')); % amp select min
        maxint = str2num(get(handles.edit_select_integral_max,'String')); % amp select max
        
        if get(handles.check_conductance,'Value')
            minint = minint./voltage;
            maxint = maxint./voltage;
            
        end
        
        int_num_points = 100;
        int_point_spacing = (xlimits(2) - xlimits(1))/int_num_points;
        int_points_vector = xlimits(1):int_point_spacing:xlimits(2);
        int_min_curve = minint./int_points_vector;
        int_max_curve = maxint./int_points_vector;
        
        int_min_plot_points = find(int_min_curve > ylimits(1) & int_min_curve < ylimits(2));
        int_max_plot_points = find(int_max_curve > ylimits(1) & int_max_curve < ylimits(2));
        
        hold on;
        plot(int_points_vector(int_min_plot_points),int_min_curve(int_min_plot_points),'-b')
        plot(int_points_vector(int_max_plot_points),int_max_curve(int_max_plot_points),'-b')
        %%%%%%%%%%%%%%%%%%%%%%
        
        % draw straight lines for dwell and amp limits in box mode
        if get(handles.radio_box_mode,'Value') == 1 % draw the box with 4 lines, if they are in the field of view
            if mindwell > xlimits(1)
                axes(handles.axes_graphs)
                hold on;
                line([mindwell mindwell],ylimits,'LineWidth',1,'Color',[0 0 1]);
            end
            if maxdwell < xlimits(2)
                axes(handles.axes_graphs)
                hold on;
                line([maxdwell maxdwell],ylimits,'LineWidth',1,'Color',[0 0 1]);
            end
            if minavg > ylimits(1)
                axes(handles.axes_graphs)
                hold on;
                line(xlimits,[minavg minavg],'LineWidth',1,'Color',[0 0 1]);
            end
            if maxavg < ylimits(2)
                axes(handles.axes_graphs)
                hold on;
                line(xlimits,[maxavg maxavg],'LineWidth',1,'Color',[0 0 1]);
            end
        end
        
        % draw the ellipse for limits
        if get(handles.radio_ellipse_mode,'Value') == 1 
            widthdwell = maxdwell - mindwell;
            heightavg = maxavg - minavg;

            axes(handles.axes_graphs)
            hold on;
            rectangle('Position',[mindwell,minavg,widthdwell,heightavg],'Curvature',[1,1]);
        end
    end
    
    % draw baseline limts if baseline plot
    if get(handles.radio_2_base,'Value') == 1 

        set(handles.text_status,'String','Drawing Baseline Limits')

        xlimits = get(handles.axes_graphs,'XLim'); % time
        ylimits = get(handles.axes_graphs,'YLim'); % baseline

        axes(handles.axes_graphs) % do it in the gui

        fun_plot_baseline(handles,'X') % update the plot

        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
        set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
        set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
        set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
        set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))

        minbase = str2num(get(handles.edit_select_baseline_min,'String')); %  select min
        maxbase = str2num(get(handles.edit_select_baseline_max,'String')); %  select max

        % draw min and max lines on the plot, if they are in the field of view
        if minbase > ylimits(1)
            axes(handles.axes_graphs)
            hold on;
            line(xlimits,[minbase minbase],'LineWidth',1,'Color',[0 0 1]);
        end
        if maxbase < ylimits(2)
            axes(handles.axes_graphs)
            hold on;
            line(xlimits,[maxbase maxbase],'LineWidth',1,'Color',[0 0 1]);
        end
    end

    % draw dwell limits on dwell hist
    if get(handles.radio_2_dwell,'Value') == 1 

        set(handles.text_status,'String','Drawing Dwell Limits')

        xlimits = get(handles.axes_graphs,'XLim'); % dwell
        ylimits = get(handles.axes_graphs,'YLim'); % count

        axes(handles.axes_graphs)
        fun_plot_dwell(handles,'X')

        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
        set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
        set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
        set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
        set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))

        mindwell = str2num(get(handles.edit_select_dwell_min,'String')); % dwell select min
        maxdwell = str2num(get(handles.edit_select_dwell_max,'String')); % dwell select max

        % draw min and max lines if in the field of view
        if mindwell > xlimits(1)
            axes(handles.axes_graphs)
            if isfield(handles,'radio2dwelllinemin')
                delete(handles.radio2dwelllinemin)
            end
            hold on;
            handles.radio2dwelllinemin = line([mindwell mindwell],ylimits,'LineWidth',1,'Color',[1 0 0]);
        end
        if maxdwell < xlimits(2)
            axes(handles.axes_graphs)
            if isfield(handles,'radio2dwelllinemax')
                delete(handles.radio2dwelllinemax)
            end
            hold on;
            handles.radio2dwelllinemax = line([maxdwell maxdwell],ylimits,'LineWidth',1,'Color',[1 0 0]);
        end

    end

    % draw amplitude limits
    if get(handles.radio_2_amplitude,'Value') == 1 

        set(handles.text_status,'String','Drawing Avg Limits')

        xlimits = get(handles.axes_graphs,'XLim'); % amp
        ylimits = get(handles.axes_graphs,'YLim'); % count

        axes(handles.axes_graphs)
        fun_plot_amplitude(handles,'X')

        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
        set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
        set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
        set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
        set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))

        minavg = str2num(get(handles.edit_select_amplitude_min,'String')); % avg select min
        maxavg = str2num(get(handles.edit_select_amplitude_max,'String')); % avg select max

        % draw min and max lines if in the field of view
        if minavg > xlimits(1)
            axes(handles.axes_graphs)
            if isfield(handles,'radio2avglinemin')
                delete(handles.radio2avglinemin)
            end
            hold on;
            handles.radio2avglinemin = line([minavg minavg],ylimits,'LineWidth',1,'Color',[1 0 0]);
        end
        if maxavg < xlimits(2)
            axes(handles.axes_graphs)
            if isfield(handles,'radio2avglinemax')
                delete(handles.radio2avglinemax)
            end
            hold on;
            handles.radio2avglinemax = line([maxavg maxavg],ylimits,'LineWidth',1,'Color',[1 0 0]);
        end
    end
    
    % draw integral (ECD) limits
    if get(handles.radio_2_integral_hist,'Value') == 1 

        set(handles.text_status,'String','Drawing Integral Limits')

        xlimits = get(handles.axes_graphs,'XLim'); % int
        ylimits = get(handles.axes_graphs,'YLim'); % counts

        axes(handles.axes_graphs)
        fun_plot_integral_hist(handles,'X')

        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
        set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
        set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
        set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
        set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))

        minint = str2num(get(handles.edit_select_integral_min,'String')); % avg select min
        maxint = str2num(get(handles.edit_select_integral_max,'String')); % avg select max

        % draw min and max lines if in the field of view
        if minint > xlimits(1)
            axes(handles.axes_graphs)
            if isfield(handles,'radio2intlinemin')
                delete(handles.radio2intlinemin)
            end
            hold on;
            handles.radio2intlinemin = line([minint minint],ylimits,'LineWidth',1,'Color',[1 0 0]);
        end
        if maxint < xlimits(2)
            axes(handles.axes_graphs)
            if isfield(handles,'radio2intlinemax')
                delete(handles.radio2intlinemax)
            end
            hold on;
            handles.radio2intlinemax = line([maxint maxint],ylimits,'LineWidth',1,'Color',[1 0 0]);
        end
    end
    
    % TODO draw limits for other plots
    
    guidata(hObject, handles); % update handles

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_mouse_x_center_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_mouse_x_center_CreateFcn(hObject, eventdata, handles)

    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_mouse_y_center_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_mouse_y_center_CreateFcn(hObject, eventdata, handles)

    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this gets the position of the mouse click in the proper units and updates
% the two text boxes, otherwise it opens the plot in a new window
function axes_graphs_ButtonDownFcn(hObject, eventdata, handles)

    set(handles.text_status,'String','Detecting Mouse Click Point')
    mouseloc = get(handles.axes_graphs,'CurrentPoint'); % button down detected
    set(handles.edit_mouse_x_center,'String',num2str(mouseloc(1,1)))
    set(handles.edit_mouse_y_center,'String',num2str(mouseloc(1,2)))

    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%   set plot min Y value   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changes y limit when y limit textbox is changed
function edit_plot2ymin_Callback(hObject, eventdata, handles)
  
    set(handles.text_status,'String','Changing Plot2 YMin')
    % get current limit
    ylimits = get(handles.axes_graphs,'YLim');
    % get new limit
    minvalue = str2double(get(handles.edit_plot2ymin,'String'));
    % only change lower limit
    set(handles.axes_graphs,'YLim',[minvalue ylimits(2)]);
    guidata(hObject, handles); % update handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_plot2ymin_CreateFcn(hObject, eventdata, handles)

    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changes y limit when y limit textbox is changed
function edit_plot2ymax_Callback(hObject, eventdata, handles)

    % set plot max Y value
    
    set(handles.text_status,'String','Changing Plot2 YMax')
    % get current limit
    ylimits = get(handles.axes_graphs,'YLim');
    % get new limit
    maxvalue = str2double(get(handles.edit_plot2ymax,'String'));
    % only change upper limit
    set(handles.axes_graphs,'YLim',[ylimits(1) maxvalue]);
    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_2_numbins_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_2_numbins_CreateFcn(hObject, eventdata, handles)

    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text_status_Callback(hObject, eventdata, handles)
% empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text_status_CreateFcn(hObject, eventdata, handles)

    if ispc
        set(hObject,'BackgroundColor','white');
    else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_use_extra_points_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot baseline
function radio_2_base_Callback(hObject, eventdata, handles)
    
    % update other radio button's states
    set(handles.text_status,'String','Plotting Baseline')
    set(handles.radio_2_scatter,'Value',0);
    set(handles.radio_2_base,'Value',1);
    set(handles.radio_2_dwell,'Value',0);
    set(handles.radio_2_amplitude,'Value',0);
    set(handles.radio_2_R_pore,'Value',0);
    set(handles.radio_2_current,'Value',0);
    set(handles.radio_2_scatter_integral,'Value',0);
    set(handles.radio_2_integral_hist,'Value',0);
    set(handles.radio_2_scatter_heat_map,'Value',0);
    set(handles.radio_2_firstextrapointmean,'Value',0);
    set(handles.radio_2_lastextrapointmean,'Value',0);

    fun_plot_baseline(handles,'X') % plot baseline
    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot baseline data
function fun_plot_baseline(handles, location) 
    % plot it here in GUI?
    if strcmp(location,'X') 
        cla(handles.axes_graphs,'reset')
        axes(handles.axes_graphs)
    end
    % open new figure?
    if strcmp(location,'n') 
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        figure()
    end

    pp3 = plot(handles.event_info_temp(:,5),'-r'); % plot baseline data

    hXLabel = xlabel('Event Number');
    hYLabel = ylabel('Baseline (nA)');
    
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
    set(pp3                     , ...
      'Color'           , [1 0 0], ...
      'LineWidth'       , 2           );

    % in GUI
    if strcmp(location,'X')

        set(pp3, 'ButtonDownFcn', {@plot2_ButtonDownFcn, handles}); % link to mouse click position function

        % get limits and update limit text boxes
        xlimits = get(handles.axes_graphs,'XLim');
        ylimits = get(handles.axes_graphs,'YLim');
        set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
        set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
        set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
        set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
    end

    if strcmp(location,'n') % use same limits
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reload dwell min & max values into the text boxes
function button_reset_dwell_Callback(hObject, eventdata, handles)

    set(handles.text_status,'String','Reseting Dwell Sort Values')

    switch get(handles.pop_time_type,'Value') 
        case 1 % FWHM_dwell
            set(handles.edit_select_dwell_min,'String',num2str(min(handles.event_info(:,4))))
            set(handles.edit_select_dwell_max,'String',num2str(max(handles.event_info(:,4))))
        case 2 % dwellold (between start and stop)
            set(handles.edit_select_dwell_min,'String',num2str(min(handles.event_more_info(:,8))))
            set(handles.edit_select_dwell_max,'String',num2str(max(handles.event_more_info(:,8))))
        otherwise
    end
    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reload baseline min & max values into the text boxes
function button_reset_baseline_Callback(hObject, eventdata, handles)

    set(handles.text_status,'String','Reseting Baseline Sort Values')
    set(handles.edit_select_baseline_min,'String',num2str(min(handles.event_info(:,5))))
    set(handles.edit_select_baseline_max,'String',num2str(max(handles.event_info(:,5))))
    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reload amplitude min & max values into the text boxes
function button_reset_amplitude_Callback(hObject, eventdata, handles)

    set(handles.text_status,'String','Reseting Avg Sort Values')
        switch get(handles.pop_amplitude_type,'Value') 
            case 1 % amplitude_minima_maxima
                set(handles.edit_select_amplitude_min,'String',num2str(min(handles.event_info(:,1))))
                set(handles.edit_select_amplitude_max,'String',num2str(max(handles.event_info(:,1))))
            case 2 % Max Amplitude
                set(handles.edit_select_amplitude_min,'String',num2str(min(handles.event_info(:,2))))
                set(handles.edit_select_amplitude_max,'String',num2str(max(handles.event_info(:,2))))
            case 3 % Integral / FWHM
                set(handles.edit_select_amplitude_min,'String',num2str(min(handles.event_more_info(:,9))))
                set(handles.edit_select_amplitude_max,'String',num2str(max(handles.event_more_info(:,9))))
            otherwise
        end
    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot amp hist
function radio_2_dwell_Callback(hObject, eventdata, handles)

    % update other radio button's states
    set(handles.text_status,'String','Plotting Dwell Hist')
    set(handles.radio_2_scatter,'Value',0);
    set(handles.radio_2_base,'Value',0);
    set(handles.radio_2_dwell,'Value',1);
    set(handles.radio_2_amplitude,'Value',0);
    set(handles.radio_2_current,'Value',0);
    set(handles.radio_2_R_pore,'Value',0);
    set(handles.radio_2_scatter_integral,'Value',0);
    set(handles.radio_2_integral_hist,'Value',0);
    set(handles.radio_2_scatter_heat_map,'Value',0);
    set(handles.radio_2_firstextrapointmean,'Value',0);
    set(handles.radio_2_lastextrapointmean,'Value',0);

    fun_plot_dwell(handles, 'X') % plot dwell in GUI
    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the dwell time histogram
function fun_plot_dwell(handles, location)
    colormap('default')
    if strcmp(location,'X') % plot here?
        cla(handles.axes_graphs,'reset')
        axes(handles.axes_graphs)
    end
    if strcmp(location,'n') % plot in new figure?
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        figure()
    end

    bins = str2num(get(handles.edit_dwellbins,'String'));
    
    switch get(handles.pop_time_type,'Value') % type of dwell data to use?
        case 1 % FWHM_dwell
            dwelldata = handles.event_info_temp(:,4);
        case 2 % dwellold (between start and stop)
            dwelldata = handles.event_more_info_temp(:,8);
        otherwise
    end

    [n,xout] = hist(dwelldata,bins); % plot it
    hist(dwelldata,bins)
    binwidth = xout(2)-xout(1);
    set(handles.edit_2_bin_width,'String',num2str(binwidth))

    hXLabel = xlabel('Time (ms)');
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

    if strcmp(location,'s') % save?
        dir_str = strcat(handles.pathname,'plots', filesep); % in the plots directory
        if(~exist(dir_str, 'dir'))
            mkdir(dir_str);
        end
        file_str = strcat(dir_str,'plot_data_dwell.dat'); % save dat file with dwell data
        dlmwrite(file_str, dwelldata, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
        
        file_str = strcat(dir_str,'plot_data_dwell_hist.dat'); % create dat file with the hist data
        dlmwrite(file_str, [xout' n'], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot amp hist
function radio_2_amplitude_Callback(hObject, eventdata, handles)

    % update other radio button's states
    set(handles.text_status,'String','Plotting Avg Hist')
    set(handles.radio_2_scatter,'Value',0);
    set(handles.radio_2_base,'Value',0);
    set(handles.radio_2_dwell,'Value',0);
    set(handles.radio_2_amplitude,'Value',1);
    set(handles.radio_2_current,'Value',0);
    set(handles.radio_2_R_pore,'Value',0);
    set(handles.radio_2_scatter_integral,'Value',0);
    set(handles.radio_2_integral_hist,'Value',0);
    set(handles.radio_2_scatter_heat_map,'Value',0);
    set(handles.radio_2_firstextrapointmean,'Value',0);
    set(handles.radio_2_lastextrapointmean,'Value',0);

    fun_plot_amplitude(handles, 'X') % plot amplitude in GUI
    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot amplitude histogram
function fun_plot_amplitude(handles, location) 
    colormap('default')
    if strcmp(location,'X') % plot in GUI?
        cla(handles.axes_graphs,'reset')
        axes(handles.axes_graphs)
    end

    if strcmp(location,'n') % plot in new figure?
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        figure()
    end

    bins = str2num(get(handles.edit_cbbins,'String'));
    voltage = str2double(get(handles.edit_voltage,'String'))/1000; % conductance case
    %%% NOTE: another copy of this code in the gaussian button func
    
    switch get(handles.pop_amplitude_type,'Value') % type of Amplitude calculation to use?
        case 1 % amplitude_minima_maxima
            if get(handles.check_conductance,'Value')
                if get(handles.check_relative,'Value')
                    amplitude = handles.event_info_temp(:,1)./handles.event_info_temp(:,5); % relative conductance case
                else
                    amplitude = handles.event_info_temp(:,1)./voltage; % conductance case
                end
            else
                if get(handles.check_relative,'Value')
                    amplitude = handles.event_info_temp(:,1)./handles.event_info_temp(:,5); % relative case
                else
                    amplitude = handles.event_info_temp(:,1); % current (nA) case
                end
                
            end        
        case 2 % Max Amplitude
            if get(handles.check_conductance,'Value')
                if get(handles.check_relative,'Value')
                    amplitude = handles.event_info_temp(:,2)./handles.event_info_temp(:,5); % relative conductance case
                else
                    amplitude = handles.event_info_temp(:,2)./voltage; % conductance case
                end
            else
                if get(handles.check_relative,'Value')
                    amplitude = handles.event_info_temp(:,2)./handles.event_info_temp(:,5); % relative case
                else
                    amplitude = handles.event_info_temp(:,2); % current (nA) case
                end
            end
        case 3 % Integral / FWHM
            if get(handles.check_conductance,'Value')
                if get(handles.check_relative,'Value')
                    amplitude = handles.event_more_info_temp(:,9)./handles.event_info_temp(:,5); % relative conductance case
                else
                    amplitude = handles.event_more_info_temp(:,9)./voltage; % conductance case
                end
            else
                if get(handles.check_relative,'Value')
                    amplitude = handles.event_more_info_temp(:,9)./handles.event_info_temp(:,5); % relative case
                else
                    amplitude = handles.event_more_info_temp(:,9); % current (nA) case
                end
            end
        otherwise
    end

    [n,xout] = hist(amplitude,bins); % plot it
    hist(amplitude,bins)
    binwidth = xout(2)-xout(1);
    set(handles.edit_2_bin_width,'String',num2str(binwidth))

    % use correct axis labels
    if get(handles.check_relative,'Value')
        hXLabel = xlabel('Relative Blockade');
    elseif get(handles.check_conductance,'Value')
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

    if strcmp(location,'s') % save?
        dir_str = strcat(handles.pathname,'plots', filesep); % in the plots dir
        if(~exist(dir_str, 'dir'))
            mkdir(dir_str);
        end
        file_str = strcat(dir_str,'plot_data_amplitude.dat'); % create dat file with the data
        dlmwrite(file_str, amplitude, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
        
        file_str = strcat(dir_str,'plot_data_amplitude_hist.dat'); % create dat file with the hist data
        dlmwrite(file_str, [xout' n'], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot current hist
function radio_2_current_Callback(hObject, eventdata, handles)
    
    % update other radio button's states
    set(handles.text_status,'String','Plotting Current Hist')
    set(handles.radio_2_scatter,'Value',0);
    set(handles.radio_2_base,'Value',0);
    set(handles.radio_2_dwell,'Value',0);
    set(handles.radio_2_amplitude,'Value',0);
    set(handles.radio_2_current,'Value',1);
    set(handles.radio_2_R_pore,'Value',0);
    set(handles.radio_2_scatter_integral,'Value',0);
    set(handles.radio_2_integral_hist,'Value',0);
    set(handles.radio_2_scatter_heat_map,'Value',0);
    set(handles.radio_2_firstextrapointmean,'Value',0);
    set(handles.radio_2_lastextrapointmean,'Value',0);

    fun_plot_current(handles,'X') % plot current hist
    guidata(hObject, handles); % update handles

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% current histogram function
function fun_plot_current(handles, location) 

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

        bins = str2num(get(handles.edit_currentbins,'String')); % number of bins

        %%%%%%%%%%%%%%%%%%%
        fun_make_new_trace(handles) % reformat trace so that baselines are zero and cutoff extrapoints if necessary
        handles=guidata(gcbo); % update local handles
        %%%%%%%%%%%%%%%%%%%

        if get(handles.check_conductance,'Value') && ~get(handles.check_relative,'Value') % convert to nS
            voltage = str2double(get(handles.edit_voltage,'String'))/1000;
            handles.newtrace = handles.newtrace./voltage;
        end


        [n,xout] = hist(handles.newtrace,bins); % make the histogram
        hist(handles.newtrace,bins)
        binwidth = xout(2)-xout(1);
        set(handles.edit_2_bin_width,'String',num2str(binwidth))

        % next four lines are needed to make the histogram show properly under
        % log scale
        ph = get(gca,'children');
        vn = get(ph,'Vertices');
        vn(:,2) = vn(:,2) + 1;
        set(ph,'Vertices',vn);

        set(gca,'YScale','log') % use log scale on Y

        if get(handles.check_relative,'Value')
            hXLabel = xlabel('Relative Blockade');
        elseif get(handles.check_conductance,'Value')
            hXLabel = xlabel('Conductance (nS)');
        else
            hXLabel = xlabel('Current Blockade (nA)');
        end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% take only the trace for selected events
function fun_make_new_trace(handles) 


    indpos = 1; % pointer to current position, initialize
    handles.extrapoints = str2num(get(handles.edit_extra_points,'String')); % get number of extra points to each side of the event
    handles.newtrace = 0;
    newstart = 0;

    number_of_selected_events = size(handles.event_info_temp,1);

    % pre allocate for speed
    newstart = zeros(number_of_selected_events);
    newend = zeros(number_of_selected_events);

    for i=1:number_of_selected_events % for number of selected events

        eventpoints = [(handles.event_info_temp(i,8)) (handles.event_info_temp(i,9))]; % get event start stop

        if get(handles.check_use_extra_points,'Value')
            eventlength = handles.event_info_temp(i,9) - handles.event_info_temp(i,8) + handles.extrapoints*2 + 1; % calc length in points
            eventtrace = handles.rawtrace_vect(eventpoints(1)-handles.extrapoints:eventpoints(2)+handles.extrapoints) - handles.event_info_temp(i,5); % get the event trace, and remove baseline
            newstart(i) = indpos + handles.extrapoints; % calculate the new start position
        else % do not use extra points
            eventlength = handles.event_info_temp(i,9) - handles.event_info_temp(i,8) + 1; % calc length in points
            eventtrace = handles.rawtrace_vect(eventpoints(1):eventpoints(2)) - handles.event_info_temp(i,5); % get the event trace, and remove baseline
            newstart(i) = indpos; % calculate the new start position
        end

        if get(handles.check_relative,'Value')
            handles.newtrace(indpos:indpos+eventlength-1) = eventtrace./handles.event_info_temp(i,5); % divide by baseline and add this trace to new trace
        else
            handles.newtrace(indpos:indpos+eventlength-1) = eventtrace; % add this trace to new trace
        end

        newend(i) = newstart(i) + handles.event_info_temp(i,9) - handles.event_info_temp(i,8); % calculate the new end position
        indpos = indpos + eventlength; % calc new start index value
    end
    guidata(gcbo, handles); % save handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot scatter
function radio_2_scatter_Callback(hObject, eventdata, handles)
    
    % update other radio button's states
    set(handles.text_status,'String','Plotting Scatter')
    set(handles.radio_2_scatter,'Value',1);
    set(handles.radio_2_base,'Value',0);
    set(handles.radio_2_dwell,'Value',0);
    set(handles.radio_2_amplitude,'Value',0);
    set(handles.radio_2_current,'Value',0);
    set(handles.radio_2_R_pore,'Value',0);
    set(handles.radio_2_scatter_integral,'Value',0);
    set(handles.radio_2_integral_hist,'Value',0);
    set(handles.radio_2_scatter_heat_map,'Value',0);
    set(handles.radio_2_firstextrapointmean,'Value',0);
    set(handles.radio_2_lastextrapointmean,'Value',0);

    fun_plot_scatter(handles,'X'); % plot the scatter plot
    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in button_plot_external.
function button_plot_external_Callback(hObject, eventdata, handles)

    %new window code here
    set(handles.text_status,'String','Generating New Figure');

    % find what type of plot
    selection = get(handles.uipanel_plot_2_type,'SelectedObject');

    % call the appropriate plotting function with the 'n' tag to indicate new
    % figure

    switch get(selection,'Tag')
        case 'radio_2_scatter'
            fun_plot_scatter(handles, 'n')
        case 'radio_2_dwell'
            fun_plot_dwell(handles, 'n')
        case 'radio_2_base'
            fun_plot_baseline(handles, 'n')
        case 'radio_2_amplitude'
            fun_plot_amplitude(handles, 'n')
        case 'radio_2_current'
            fun_plot_current(handles, 'n')
        case 'radio_2_R_pore'
            fun_plot_R_pore(handles, 'n')
        case 'radio_2_scatter_integral'
            fun_plot_scatter_integral(handles, 'n')
        case 'radio_2_integral_hist'
            fun_plot_integral_hist(handles, 'n')
        case 'radio_2_scatter_heat_map'
            fun_plot_scatter_heat(handles, 'n')
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_multiply_by_ten_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dump variables to Matlab base workspace
function button_dump_variables_Callback(hObject, eventdata, handles)
    
    % TODO update, missing newer variables?
    if isfield(handles,'gui_events_file_columns')
        assignin('base', 'gui_events_file_columns', handles.gui_events_file_columns)
    end
    if isfield(handles,'extrapoints')
        assignin('base', 'extrapoints', handles.extrapoints)
    end
    if isfield(handles,'timestep')
        assignin('base', 'timestep', handles.timestep)
    end
    if isfield(handles,'event_info')
        assignin('base', 'event_info', handles.event_info)
    end
    if isfield(handles,'event_info_temp')
        assignin('base', 'event_info_temp', handles.event_info_temp)
    end
    if isfield(handles,'rawtrace_vect')
        assignin('base', 'rawtrace_vect', handles.rawtrace_vect)
    end
    if isfield(handles,'event_more_info')
        assignin('base', 'event_more_info', handles.event_more_info)
    end
    if isfield(handles,'event_more_info_temp')
        assignin('base', 'event_more_info_temp', handles.event_more_info_temp)
    end
    if isfield(handles,'rawunfilteredtrace_vect')
        assignin('base', 'rawunfilteredtrace_vect', handles.rawunfilteredtrace_vect)
    end
    if isfield(handles,'newtrace')
        assignin('base', 'newtrace', handles.newtrace)
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Opens the Gaussian fitting GUI and passes the necessary data
function button_gaussian_fit_Callback(hObject, eventdata, handles)

    %
    % Important!!!
    %
    % there is a copy of this code in GUI_localstructures
    % in the function button_fit_gaussian_Callback
    % If you change it here, change it there too.

    HandleMainGUI=getappdata(0,'HandleMainGUI'); % handle used to pass data to new GUI

    switch get(handles.pop_time_type,'Value') 
        case 1 % FWHM_dwell
            setappdata(HandleMainGUI,'dwellvector',handles.event_info_temp(:,4)); % pass dwell data to sub GUI
        case 2 % dwellold (between start and stop)
            setappdata(HandleMainGUI,'dwellvector',handles.event_more_info_temp(:,8)); % pass dwell data to sub GUI
        otherwise
    end


    if isfield(handles,'rawtrace_vect')
        fun_make_new_trace(handles) % make new trace with only selected events
        handles=guidata(gcbo); % grab new trace
        setappdata(HandleMainGUI,'newtrace',handles.newtrace); % pass it down
    end

    setappdata(HandleMainGUI,'dwellbins',str2num(get(handles.edit_dwellbins,'String'))); % pass number of dwell bins
    setappdata(HandleMainGUI,'cbbins',str2num(get(handles.edit_cbbins,'String'))); % pass number of amplitude bins
    setappdata(HandleMainGUI,'currentbins',str2num(get(handles.edit_currentbins,'String'))); % pass number of current bins
    
    voltage = str2double(get(handles.edit_voltage,'String'))/1000; % get voltage
    amplitude = fun_get_amplitude(handles);

    setappdata(HandleMainGUI,'amplitude',amplitude); % pass it down

    % tell the GUI what method you use for amplitude
    if get(handles.check_relative,'Value')
        check_conductance = 2;
    elseif get(handles.check_conductance,'Value')
        check_conductance = 1;
    else
        check_conductance = 0;
    end
    
    setappdata(HandleMainGUI,'check_conductance',check_conductance); % use conductance or current or relative
    voltage = abs(str2double(get(handles.edit_voltage,'String'))/1000); 
    
    setappdata(HandleMainGUI,'pathname',handles.pathname);
    GUI_gaussian % open gaussian fitting GUI

    guidata(hObject, handles); % update handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this plots all of the different types of plots available and saves them
% in .fig, .eps, and .png format in the \plots directory
function button_mass_export_Callback(hObject, ~, handles)

    set(handles.text_status,'String','Mass exporting all plots to events file directory....')

    dir_str = strcat(handles.pathname,'plots', filesep); % in this sub directory
    if(~exist(dir_str, 'dir'))
        mkdir(dir_str);
    end

    % scatter plot
    newFig = figure;
    fun_plot_scatter(handles, 's') % data file saved inside this function
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_scatter.eps')) 
    if exist('export_fig') == 2
        set(newFig, 'Color', 'w');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter_exp_fig.eps'), '-cmyk');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter_exp_fig.png'), '-cmyk');
    end
    close(newFig)

    % scatter with hist plot
    newFig = figure;
    fun_plot_scatter(handles, 'h') % data file saved inside this function
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter-hist.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter-hist_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter-hist_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_scatter-hist.eps')) 
    if exist('export_fig') == 2
        set(newFig, 'Color', 'w');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter-hist_exp_fig.eps'), '-cmyk');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter-hist_exp_fig.png'), '-cmyk');
    end
    close(newFig)

    % scatter heat
    newFig = figure;
    fun_plot_scatter_heat(handles, 's') % data file saved inside this function
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter_heat.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter_heat_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter_heat_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_scatter_heat.eps')) 
    if exist('export_fig') == 2
        set(newFig, 'Color', 'w');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter_heat_exp_fig.eps'), '-cmyk');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter_heat_exp_fig.png'), '-cmyk');
    end
    close(newFig)

    
    % dwell hist plot
    newFig = figure;
    fun_plot_dwell(handles, 's') % data file saved inside this function
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_dwell_hist.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_dwell_hist_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_dwell_hist_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_dwell_hist.eps')) 
    if exist('export_fig') == 2
        set(newFig, 'Color', 'w');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_dwell_hist_exp_fig.eps'), '-cmyk');
        export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_dwell_hist_exp_fig.png'), '-cmyk');
    end
    close(newFig)

    % baseline plot
    newFig = figure;
    fun_plot_baseline(handles, 's') 
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_baseline.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_baseline_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_baseline_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_baseline.eps')) 
    close(newFig)

    % amplitude histogram plot
    newFig = figure;
    fun_plot_amplitude(handles, 's') % data file saved inside this function
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_condblockade_hist.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_condblockade_hist_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_condblockade_hist_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_condblockade.eps')) 
    close(newFig)

    % current hist plot
    if isfield(handles,'rawtrace_vect') % if rawtrace is loaded
        newFig = figure;
        fun_plot_current(handles, 's') % data file saved inside this function
        saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_current_hist.fig'), 'fig'); % save as .fig
        saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_current_hist_large.png'), 'png');
        set(newFig, 'PaperPositionMode', 'auto');
        saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_current_hist_small.png'), 'png');
        print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_current_hist.eps')) 
        if exist('export_fig') == 2
            set(newFig, 'Color', 'w');
            export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_current_hist_exp_fig.eps'), '-cmyk');
            export_fig(newFig, strcat(handles.pathname,'plots', filesep,'plot_current_hist_exp_fig.png'), '-cmyk');
        end
        close(newFig)
    end

    % Pore resistence as a function of time plot 
    % note this is not changed by any applied selection, it is always the same
    if isfield(handles,'event_more_info') % if extra info exits, legacy check
        newFig = figure;
        fun_plot_R_pore(handles, 's')
        saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_Rpore.fig'), 'fig'); % save as .fig
        saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_Rpore_large.png'), 'png');
        set(newFig, 'PaperPositionMode', 'auto');
        saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_Rpore_small.png'), 'png');
        print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_Rpore.eps')) 
        close(newFig)
    end

    %%%%%%%%%%%%%%%%%%%
    % make a giant plot with 4 subplots of scatter, amplitude, dwell, and
    % either current or baseline
    %%%%%%%%%%%%%%%%%%%
    newFig = figure;

    subplot(2,2,1)
    fun_plot_scatter(handles, 's')

    subplot(2,2,2)
    fun_plot_amplitude(handles, 's')

    subplot(2,2,3)
    fun_plot_dwell(handles, 's')

    subplot(2,2,4)
    if isfield(handles,'rawtrace_vect')
        fun_plot_current(handles, 's')
    else
        fun_plot_baseline(handles, 's')
    end

    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_main_all.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_main_all_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_main_all_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_main_all.eps'))
    close(newFig)
    %%%%%%%%%%%%%%%%% end massive plot

    % plot scatter int
    newFig = figure;
    fun_plot_scatter_integral(handles, 's') % data file saved inside this function
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter_integral.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter_integral_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_scatter_integral_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_scatter_integral.eps'))
    close(newFig)

    % plot int hist
    newFig = figure;
    fun_plot_integral_hist(handles, 's')
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_integral_hist.fig'), 'fig'); % save as .fig
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_integral_hist_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_integral_hist_small.png'), 'png');
    print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_integral_hist.eps'))
    close(newFig)

    % all done here

    set(handles.text_status,'String','Mass exporting all plots to events file directory....done')

    guidata(hObject, handles); % update handles



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_amplitude_type_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_amplitude_type_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_conductance_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_open_event_plot_Callback(hObject, eventdata, handles)

    fun_plot_event('n') % open it in new figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_show_baseline_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_voltage_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_voltage_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_time_type_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_time_type_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_overlay_unfiltered_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_overlay_unfiltered_CreateFcn(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in button_ghost.
function button_ghost_Callback(hObject, eventdata, handles)

    % prepare data for opening the ghosting GUI
    HandleMainGUI=getappdata(0,'HandleMainGUI'); % make a master handle

    eventnumber = str2double(get(handles.text_event_selected,'String')); % find the event to ghost

    eventpoints = [(handles.event_info_temp(eventnumber,8)) (handles.event_info_temp(eventnumber,9))]; % event start stop

    handles.extrapoints = str2num(get(handles.edit_extra_points,'String')); % # of extra points

    eventtrace = handles.rawtrace_vect(eventpoints(1)-handles.extrapoints:eventpoints(2)+handles.extrapoints)-handles.event_info_temp(eventnumber,5); % build the trace
    % build the time vector
    timetemp = transpose(0:(handles.timestep*1000):handles.timestep*(eventpoints(2)-eventpoints(1)+2*handles.extrapoints)*1000);

    if isfield(handles,'rawunfilteredtrace_vect') % ghost should not be enabled if unfiltered is not loaded but try to catch it anyway
        eventtraceunfilt = handles.rawunfilteredtrace_vect(eventpoints(1)-handles.extrapoints:eventpoints(2)+handles.extrapoints)-handles.event_info_temp(eventnumber,5);
    end

    eventlengthinpoints = eventpoints(2) - eventpoints(1) +1;

    % pass all data to ghost GUI
    setappdata(HandleMainGUI,'eventtrace',eventtrace);
    setappdata(HandleMainGUI,'timetemp',timetemp);
    setappdata(HandleMainGUI,'eventtraceunfilt',eventtraceunfilt);
    setappdata(HandleMainGUI,'extrapoints',handles.extrapoints);
    setappdata(HandleMainGUI,'eventlengthinpoints',eventlengthinpoints);
    setappdata(HandleMainGUI,'timestep',handles.timestep);
    setappdata(HandleMainGUI,'fsample',handles.fsample);

    GUI_ghost % Who you gonna call?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_ghost_CreateFcn(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in button_event_rate_mgmt.
function button_event_rate_mgmt_Callback(hObject, eventdata, handles)
    % get ready to run the event rate GUI

    % send localtracetime, totaltime, numevents, extra_info_temp

    % this was useful in learning to pass data to another GUI
    %http://www.mathworks.com/matlabcentral/answers/338-how-to-pass-data-from-o
    %ne-gui-to-another

    HandleMainGUI=getappdata(0,'HandleMainGUI'); % make handle

    % send data
    setappdata(HandleMainGUI,'localtracetime',handles.localtracetime);
    switch get(handles.pop_time_type,'Value') 
        case 1 % FWHM_dwell
            setappdata(HandleMainGUI,'dwellvector',handles.event_info_temp(:,4)); % pass dwell data to sub GUI
        case 2 % dwellold (between start and stop)
            setappdata(HandleMainGUI,'dwellvector',handles.event_more_info_temp(:,8)); % pass dwell data to sub GUI
        otherwise
    end
    setappdata(HandleMainGUI,'pathname',handles.pathname);
    setappdata(HandleMainGUI,'totaltimetrace',handles.totaltimetrace);
    setappdata(HandleMainGUI,'totaleventsfound',handles.totaleventsfound);
    setappdata(HandleMainGUI,'event_more_info_temp',handles.event_more_info_temp);

    GUI_eventrate % Run GUI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GUI_events_CloseRequestFcn(hObject, eventdata, handles)

    % clean up shared data
    HandleMainGUI=getappdata(0,'HandleMainGUI');

    % Remember that your GUIs might try to getappdata that doesn't exist, you should first test if it does exist
    if (isappdata(0,'HandleMainGUI') && isappdata(HandleMainGUI,'MySharedData'))
        rmappdata(HandleMainGUI,'MySharedData') % do rmappdata for all data shared 
    else
        %do something else, maybe loading default values into those variables
    end

    % closes the figure
    delete(hObject);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_event_rate_mgmt_CreateFcn(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this plots Rpore in nA or Cpore in nS, button 
function radio_2_R_pore_Callback(hObject, eventdata, handles) % this plots Rpore in nA or Cpore in nS

    % update other radio button's states
    if isfield(handles,'event_more_info') % legacy catch
        set(handles.text_status,'String','Plotting Rpore')
        set(handles.radio_2_scatter,'Value',0);
        set(handles.radio_2_base,'Value',0);
        set(handles.radio_2_dwell,'Value',0);
        set(handles.radio_2_amplitude,'Value',0);
        set(handles.radio_2_current,'Value',0);
        set(handles.radio_2_R_pore,'Value',1);
        set(handles.radio_2_scatter_integral,'Value',0);
        set(handles.radio_2_integral_hist,'Value',0);
        set(handles.radio_2_scatter_heat_map,'Value',0);
        set(handles.radio_2_firstextrapointmean,'Value',0);
        set(handles.radio_2_lastextrapointmean,'Value',0);

        fun_plot_R_pore(handles, 'X') % call plotting function
    end
    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this plots Rpore in nA or Cpore in nS
function fun_plot_R_pore(handles, location) 
    
    
    
    if strcmp(location,'X') % plot it here
        cla(handles.axes_graphs,'reset')
        axes(handles.axes_graphs)
    end

    if strcmp(location,'n') % plot it in a new figure
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        figure()
    end

    % this is not elegant code
    % the problem is that R pore is only calculated once per trace analyzed but
    % passed as a variable for all events in that trace
    % but we only need one R value per trace, so this extracts just one value
    % from each trace, and ignores the duplicates

    
    averagecounter = 1;
    alleventcount = 1;

    filenumbers = unique(handles.event_more_info_temp(:,1), 'rows'); % total # of trace files detected

    for filecounter=1:length(filenumbers) % for each trace
        currentfile = filenumbers(filecounter);

        eventindexes = find(handles.event_more_info_temp(:,1) == currentfile); % get Rpore value

        if ~isempty(eventindexes) % catch case of no events

            fileresistance(filecounter) = handles.event_more_info_temp(eventindexes(1),2); % add the Rpore value
            
            fileresistance_from_baseline(filecounter) = handles.voltage/mean(handles.event_info_temp(eventindexes,5)); % average of all event's baselines in that file
            
            filetime(filecounter) = handles.event_more_info_temp(eventindexes(1),4); % time of experiment for that value

        end
    end

    pp3 = plot(filetime,fileresistance,'-r'); % plot Rpore vs experimental time
    hold on
    plot(filetime,fileresistance_from_baseline,'-b');
    hXLabel = xlabel('Time (s)');
    hYLabel = ylabel('Resistance (Mohm)');
    
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
    set(pp3                     , ...
      'Color'           , [1 0 0], ...
      'LineWidth'       , 2           );

    if strcmp(location,'X')
        set(pp3, 'ButtonDownFcn', {@plot2_ButtonDownFcn, handles}); % link to mouse click position function
        % get limits and update limit text boxes
        xlimits = get(handles.axes_graphs,'XLim');
        ylimits = get(handles.axes_graphs,'YLim');
        set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
        set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
        set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
        set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
    end
    if strcmp(location,'n') % set limits for new figure
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
    end
    if strcmp(location,'s') % save
        dir_str = strcat(handles.pathname,'plots', filesep); % in the plots directory
        if(~exist(dir_str, 'dir'))
            mkdir(dir_str);
        end
        file_str = strcat(dir_str,'plot_data_Rpore.dat'); % save dat file with Rpore data
        dataoutput = [filetime', fileresistance', fileresistance_from_baseline'];
        dlmwrite(file_str, dataoutput, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_load_rawtrace_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_event_number_min_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_event_number_min_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_event_number_max_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_event_number_max_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reload event # min & max values into the text boxes
function button_reset_event_number_Callback(hObject, eventdata, handles)  
    
    set(handles.text_status,'String','Reseting Event Number Sort Values')
    set(handles.edit_select_event_number_min,'String','1')
    set(handles.edit_select_event_number_max,'String',num2str(length(handles.event_info(:,5))))
    guidata(hObject, handles); % update handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_integral_min_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_integral_min_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_integral_max_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_integral_max_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_reset_integral_Callback(hObject, eventdata, handles)

    % reload int min & max values into the text boxes
    
    set(handles.text_status,'String','Reseting Int Sort Values')
    set(handles.edit_select_integral_min,'String',num2str(min(handles.event_info(:,3))))
    set(handles.edit_select_integral_max,'String',num2str(max(handles.event_info(:,3))))
    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot scatter integral hist button
function radio_2_scatter_integral_Callback(hObject, eventdata, handles)
    
    % update other radio button's states
    set(handles.text_status,'String','Plotting Scatter Int')
    set(handles.radio_2_scatter,'Value',0);
    set(handles.radio_2_base,'Value',0);
    set(handles.radio_2_dwell,'Value',0);
    set(handles.radio_2_amplitude,'Value',0);
    set(handles.radio_2_current,'Value',0);
    set(handles.radio_2_R_pore,'Value',0);
    set(handles.radio_2_scatter_integral,'Value',1);
    set(handles.radio_2_integral_hist,'Value',0);
    set(handles.radio_2_scatter_heat_map,'Value',0);
    set(handles.radio_2_firstextrapointmean,'Value',0);
    set(handles.radio_2_lastextrapointmean,'Value',0);

    fun_plot_scatter_integral(handles,'X'); % call scatter int plot function
    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this plots scatter integral
function fun_plot_scatter_integral(handles,location) 

    if strcmp(location,'X') % plot it here
        cla(handles.axes_graphs,'reset')
        axes(handles.axes_graphs)
    end

    if strcmp(location,'n') % plot in new figure
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        figure()
    end

    voltage = str2double(get(handles.edit_voltage,'String'))/1000; % get voltage

    integrals = handles.event_info_temp(:,3); % get int data

    % get correct dwell data
    switch get(handles.pop_time_type,'Value') 
        case 1 % FWHM_dwell
            dwelldata = handles.event_info_temp(:,4);
        case 2 % dwellold (between start and stop)
            dwelldata = handles.event_more_info_temp(:,8);
        otherwise
    end

    pp3 = scatter(dwelldata,integrals,'or'); % plot it

    hXLabel = xlabel('Time (ms)');
    hYLabel = ylabel('Integral (nA*ms)');

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
  
    if get(handles.check_dwell_X_log,'Value')
        set( gca                       , ...
        'XScale'   , 'log' );
    end
    if strcmp(location,'n') % set correct limits
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
    end

    if strcmp(location,'X')
        set(pp3, 'ButtonDownFcn', {@plot2_ButtonDownFcn, handles}); % link to mouse click position function
        % get limits and update limit text boxes
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        set(handles.edit_plot2xmin,'String',num2str(xlimits(1)))
        set(handles.edit_plot2xmax,'String',num2str(xlimits(2)))
        set(handles.edit_plot2ymin,'String',num2str(ylimits(1)))
        set(handles.edit_plot2ymax,'String',num2str(ylimits(2)))
    end        

    if strcmp(location,'s') % save
        dir_str = strcat(handles.pathname,'plots', filesep); % in the plots directory
        if(~exist(dir_str, 'dir'))
            mkdir(dir_str);
        end
        file_str = strcat(dir_str,'plot_data_scatter_integral.dat'); % make a .dat file with the data needed
        datatooutput = [dwelldata,integrals];
        dlmwrite(file_str, datatooutput, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot ECD (integral) hist button
function radio_2_integral_hist_Callback(hObject, eventdata, handles)

    % update other radio button's states
    set(handles.text_status,'String','Plotting Int Hist')
    set(handles.radio_2_scatter,'Value',0);
    set(handles.radio_2_base,'Value',0);
    set(handles.radio_2_dwell,'Value',0);
    set(handles.radio_2_amplitude,'Value',0);
    set(handles.radio_2_current,'Value',0);
    set(handles.radio_2_R_pore,'Value',0);
    set(handles.radio_2_scatter_integral,'Value',0);
    set(handles.radio_2_integral_hist,'Value',1);
    set(handles.radio_2_scatter_heat_map,'Value',0);
    set(handles.radio_2_firstextrapointmean,'Value',0);
    set(handles.radio_2_lastextrapointmean,'Value',0);
    
    fun_plot_integral_hist(handles,'X'); % call int hist plotting  function
    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots integral histogram
function fun_plot_integral_hist(handles, location) 

    
    colormap('default')
    
    if strcmp(location,'X') % plot it here in the GUI
        cla(handles.axes_graphs,'reset')
        axes(handles.axes_graphs)
    end

    if strcmp(location,'n') % plot it in a new figure
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        figure()
    end

    bins = str2num(get(handles.edit_intbins,'String')); % number of bins


    [n,xout] = hist(handles.event_info_temp(:,3),bins); % plot it
    hist(handles.event_info_temp(:,3),bins)
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

    if strcmp(location,'n') % set correct limits
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
    end
    
    if strcmp(location,'s')
        dir_str = strcat(handles.pathname,'plots', filesep); % in the plots dir
        if(~exist(dir_str, 'dir'))
            mkdir(dir_str);
        end
        file_str = strcat(dir_str,'plot_data_integral.dat'); % creat dat file with the data
        dlmwrite(file_str, handles.event_info_temp(:,3), 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
        file_str = strcat(dir_str,'plot_data_integral_hist.dat'); % creat dat file with the data
        dlmwrite(file_str, [xout' n'], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in button_save_salt.
function button_save_salt_Callback(hObject, eventdata, handles)

    % call salt save function
    fun_save_salt(handles,1,pwd)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this saves salt files with all current parameter values
function fun_save_salt(handles,askgui,pathvar) 

    % askgui will either prompt user for a filename or otherwise use the
    % default name
    if askgui % get file location
        [filename, pathname] = uiputfile(...
         {'*.salt';'*.*'},...
         'Save Salt as');
    else % use default name if the salt is clean.salt
        pathname = pathvar;
        filename = 'clean.salt';
    end
    
    dir_str = strcat(pathname, filesep);
    file_str = strcat(dir_str,filename);
    fid = fopen(file_str, 'w');
    if fid == -1
        % could not open file for writing
        set(handles.text_status,'String','Could not open SALT file for writing.')
    else % write all the parameters
        handles=guidata(gcbo);
        fprintf(fid, '%s\r\n', '%% Events Salt File');
        fprintf(fid, '%s\r\n', '%% Dwell Var');
        fprintf(fid, '%s\r\n', num2str(get(handles.pop_time_type,'Value')));
        fprintf(fid, '%s\r\n', '%% Amplitude Var');
        fprintf(fid, '%s\r\n', num2str(get(handles.pop_amplitude_type,'Value')));
        fprintf(fid, '%s\r\n', '%% Voltage');
        fprintf(fid, '%s\r\n', get(handles.edit_voltage,'String'));
        fprintf(fid, '%s\r\n', '%% Extra Points');
        fprintf(fid, '%s\r\n', get(handles.edit_extra_points,'String'));
        fprintf(fid, '%s\r\n', '%% Bins');
        fprintf(fid, '%s\r\n', get(handles.edit_2_numbins,'String'));
        fprintf(fid, '%s\r\n', '%% Box-1 or Ellipse-0');
        fprintf(fid, '%s\r\n', num2str(get(handles.radio_box_mode,'Value')));
        switch get(handles.pop_time_type,'Value') 
            case 1 % FWHM_dwell
                fprintf(fid, '%s\r\n', '%% Dwell Min');
                %fprintf(fid, '%s\r\n', num2str(min(handles.event_info_temp(:,4))));
                fprintf(fid, '%s\r\n', get(handles.edit_select_dwell_min,'String'));
                fprintf(fid, '%s\r\n', '%% Dwell Max');
                %fprintf(fid, '%s\r\n', num2str(max(handles.event_info_temp(:,4))));
                fprintf(fid, '%s\r\n', get(handles.edit_select_dwell_max,'String'));
            case 2 % dwellold (between start and stop)
                fprintf(fid, '%s\r\n', '%% Dwell Min');
                %fprintf(fid, '%s\r\n', num2str(min(handles.event_more_info_temp(:,8))));
                fprintf(fid, '%s\r\n', get(handles.edit_select_dwell_min,'String'));
                fprintf(fid, '%s\r\n', '%% Dwell Max');
                %fprintf(fid, '%s\r\n', num2str(max(handles.event_more_info_temp(:,8))));
                fprintf(fid, '%s\r\n', get(handles.edit_select_dwell_max,'String'));
            otherwise
        end
        switch get(handles.pop_amplitude_type,'Value')
            case 1 % amplitude_minima_maxima
                fprintf(fid, '%s\r\n', '%% Amp Min');
                %fprintf(fid, '%s\r\n', num2str(min(handles.event_info_temp(:,1))));
                fprintf(fid, '%s\r\n', get(handles.edit_select_amplitude_min,'String'));
                fprintf(fid, '%s\r\n', '%% Amp Max');
                %fprintf(fid, '%s\r\n', num2str(max(handles.event_info_temp(:,1))));
                fprintf(fid, '%s\r\n', get(handles.edit_select_amplitude_max,'String'));
            case 2 % Max Amplitude
                fprintf(fid, '%s\r\n', '%% Amp Min');
                %fprintf(fid, '%s\r\n', num2str(min(handles.event_info_temp(:,2))));
                fprintf(fid, '%s\r\n', get(handles.edit_select_amplitude_min,'String'));
                fprintf(fid, '%s\r\n', '%% Amp Max');
                %fprintf(fid, '%s\r\n', num2str(max(handles.event_info_temp(:,2))));
                fprintf(fid, '%s\r\n', get(handles.edit_select_amplitude_max,'String'));
            case 3 % Integral / FWHM
                fprintf(fid, '%s\r\n', '%% Amp Min');
                %fprintf(fid, '%s\r\n', num2str(min(handles.event_more_info_temp(:,9))));
                fprintf(fid, '%s\r\n', get(handles.edit_select_amplitude_min,'String'));
                fprintf(fid, '%s\r\n', '%% Amp Max');
                %fprintf(fid, '%s\r\n', num2str(max(handles.event_more_info_temp(:,9))));
                fprintf(fid, '%s\r\n', get(handles.edit_select_amplitude_max,'String'));
            otherwise
        end
        fprintf(fid, '%s\r\n', '%% Base Min');
        %fprintf(fid, '%s\r\n', num2str(min(handles.event_info_temp(:,5))));
        fprintf(fid, '%s\r\n', get(handles.edit_select_baseline_min,'String'));
        fprintf(fid, '%s\r\n', '%% Base Max');
        %fprintf(fid, '%s\r\n', num2str(max(handles.event_info_temp(:,5))));
        fprintf(fid, '%s\r\n', get(handles.edit_select_baseline_max,'String'));
        fprintf(fid, '%s\r\n', '%% Event Num Min');
        fprintf(fid, '%s\r\n', get(handles.edit_select_event_number_min,'String'));
        fprintf(fid, '%s\r\n', '%% Event Num Max');
        fprintf(fid, '%s\r\n', get(handles.edit_select_event_number_max,'String'));
        fprintf(fid, '%s\r\n', '%% Int Min');
        %fprintf(fid, '%s\r\n', num2str(min(handles.event_info_temp(:,3))));
        fprintf(fid, '%s\r\n', get(handles.edit_select_integral_min,'String'));
        fprintf(fid, '%s\r\n', '%% Int Max');
        %fprintf(fid, '%s\r\n', num2str(max(handles.event_info_temp(:,3))));
        fprintf(fid, '%s\r\n', get(handles.edit_select_integral_max,'String'));
        fprintf(fid, '%s\r\n', '%% unused');
        fprintf(fid, '%s\r\n', num2str(0));
        fprintf(fid, '%s\r\n', '%% unused');
        fprintf(fid, '%s\r\n', num2str(0));
        fprintf(fid, '%s\r\n', '%% unused');
        fprintf(fid, '%s\r\n', num2str(0));
        fprintf(fid, '%s\r\n', '%% unused');
        fprintf(fid, '%s\r\n', num2str(0));
        fprintf(fid, '%s\r\n', '%% First Extra Points Mean Max');
        fprintf(fid, '%s\r\n', get(handles.edit_select_firstextrapointmeanmax,'String'));
        fprintf(fid, '%s\r\n', '%% Last Extra Points Mean Max');
        fprintf(fid, '%s\r\n', get(handles.edit_select_lastextrapointmeanmax,'String'));
        fprintf(fid, '%s\r\n', '%% END SALT FILE %%');
    end % done
    fclose(fid); % close file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_load_salt_Callback(hObject, eventdata, handles)

    % load a SALT file
    fun_load_salt(handles,1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load a SALT file, update all of the text boxes
function fun_load_salt(handles,askgui)

    if askgui
        [FileName,PathName] = uigetfile('*.salt','Select the MATLAB code file'); % get SALT location from user
    else
        PathName = handles.pathname;
        FileName = 'clean.salt';
    end
    allname = strcat(PathName, filesep, FileName); % full path and filename

    if isequal(exist(allname),2) % check it exists
        fid = fopen(allname); % open it
        if fid == -1
            disp(strcat(FileName,' SALT file could not be loaded!!'))
        else % read in parameters and set appropriate text boxes
            line = fgetl(fid); %% Events Salt File
            line = fgetl(fid); %% Dwell Var
            tempval = str2double(fgetl(fid)); %% 
            set(handles.pop_time_type,'Value',tempval);     
            line = fgetl(fid); %% Amplitude Var
            tempval = str2double(fgetl(fid)); %% 
            set(handles.pop_amplitude_type,'Value',tempval); 
            line = fgetl(fid); %% Voltage
            tempval = str2double(fgetl(fid)); %% 
            set(handles.edit_voltage,'String',num2str(tempval)); 
            line = fgetl(fid); %% Extra Points
            tempval = str2double(fgetl(fid)); %% 
            set(handles.edit_extra_points,'String',num2str(tempval)); 
            line = fgetl(fid); %% Bins
            tempval = str2double(fgetl(fid)); %% 
            set(handles.edit_2_numbins,'String',num2str(tempval));
            line = fgetl(fid); %% Box-1 or Ellipse-0
            tempval = str2double(fgetl(fid)); %% 
            set(handles.radio_box_mode,'Value',tempval);
            line = fgetl(fid); %% Dwell Min
            tempval = str2double(fgetl(fid)); %% 
            set(handles.edit_select_dwell_min,'String',num2str(tempval));
            line = fgetl(fid); %% Dwell Max
            tempval = str2double(fgetl(fid)); %% 
            set(handles.edit_select_dwell_max,'String',num2str(tempval));
            line = fgetl(fid); %% Amp Min
            tempval = str2double(fgetl(fid)); %% 
            set(handles.edit_select_amplitude_min,'String',num2str(tempval));
            line = fgetl(fid); %% Amp Max
            tempval = str2double(fgetl(fid)); %% 
            set(handles.edit_select_amplitude_max,'String',num2str(tempval));
            line = fgetl(fid); %% Base Min
            tempval = str2double(fgetl(fid)); %% 
            set(handles.edit_select_baseline_min,'String',num2str(tempval));
            line = fgetl(fid); %% Base Max
            tempval = str2double(fgetl(fid)); %% 
            set(handles.edit_select_baseline_max,'String',num2str(tempval));
            line = fgetl(fid); %% Event Num Min
            tempval = str2double(fgetl(fid)); %% 
            set(handles.edit_select_event_number_min,'String',num2str(tempval));
            line = fgetl(fid); %% Event Num Max
            tempval = str2double(fgetl(fid)); %% 
            set(handles.edit_select_event_number_max,'String',num2str(tempval));
            line = fgetl(fid); %% Int Min
            tempval = str2double(fgetl(fid)); %% 
            set(handles.edit_select_integral_min,'String',num2str(tempval));
            line = fgetl(fid); %% Int Max
            tempval = str2double(fgetl(fid)); %% 
            set(handles.edit_select_integral_max,'String',num2str(tempval));
            line = fgetl(fid); %% unused
            line = fgetl(fid); %% #
            line = fgetl(fid); %% unused
            line = fgetl(fid); %% #
            line = fgetl(fid); %% unused
            line = fgetl(fid); %% #
            line = fgetl(fid); %% unused
            line = fgetl(fid); %% #
            line = fgetl(fid); %% EOF or First Extra Points Mean Max
            if(~strcmp(line,'%% END SALT FILE %%'))
                tempval = str2double(fgetl(fid)); %% 
                set(handles.edit_select_firstextrapointmeanmax,'String',num2str(tempval));
                line = fgetl(fid); %% Last Extra Points Mean Max
                tempval = str2double(fgetl(fid)); %% 
                set(handles.edit_select_lastextrapointmeanmax,'String',num2str(tempval));
            end
        end
        fclose(fid);
    else % something went wrong
        set(handles.text_status,'String','No Salt file found')
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_lock_baseline_variation_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_baseline_variation_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_baseline_variation_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
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
function button_save_clean_Callback(hObject, eventdata, handles)

    fun_save_salt(handles,0,handles.pathname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_load_clean_Callback(hObject, eventdata, handles)

    fun_load_salt(handles,0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_save_event_Callback(hObject, eventdata, handles)

    fun_plot_event('s') % open it in new figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_dwell_X_log_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_fold_min_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_fold_min_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_fold_max_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_fold_max_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_reset_fold_count_Callback(hObject, eventdata, handles)

    set(handles.text_status,'String','Reseting Event Current Sort Values')
    set(handles.edit_select_fold_min,'String',num2str(min(handles.events_minmax(:,1))))
    set(handles.edit_select_fold_max,'String',num2str(max(handles.events_minmax(:,2))))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_use_fold_counts_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uipanel_plot_2_type_ButtonDownFcn(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% radiobutton to plot heat map
function radio_2_scatter_heat_map_Callback(hObject, eventdata, handles)

    % update other radio button's states
    set(handles.text_status,'String','Plotting Int Hist')
    set(handles.radio_2_scatter,'Value',0);
    set(handles.radio_2_base,'Value',0);
    set(handles.radio_2_dwell,'Value',0);
    set(handles.radio_2_amplitude,'Value',0);
    set(handles.radio_2_current,'Value',0);
    set(handles.radio_2_R_pore,'Value',0);
    set(handles.radio_2_scatter_integral,'Value',0);
    set(handles.radio_2_integral_hist,'Value',0);
    set(handles.radio_2_scatter_heat_map,'Value',1);
    set(handles.radio_2_firstextrapointmean,'Value',0);
    set(handles.radio_2_lastextrapointmean,'Value',0);

    fun_plot_scatter_heat(handles,'X'); % call heat map plotting  function
    guidata(hObject, handles); % update handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_dwellbins_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_dwellbins_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_cbbins_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_cbbins_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_currentbins_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_currentbins_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_intbins_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_intbins_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_heatbins_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_heatbins_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit15_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this updates the detection level text box
function fun_update_detection_level(handles)

    voltage = abs(str2double(get(handles.edit_voltage,'String'))/1000); % get voltage
    detection_levels_nS = abs(abs(handles.event_info_temp(:,5))-abs(handles.event_info_temp(:,6)))./voltage;
    % abs(abs(baseline)-abs(detection level))/voltage
    % in nS
    detection_levels_nS_avg = mean(detection_levels_nS);
    detection_levels_nS_std = std(detection_levels_nS);

    set(handles.edit_avg_detect_level,'String',num2str(detection_levels_nS_avg))
    set(handles.edit_avg_detect_level_std,'String',num2str(detection_levels_nS_std))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_num_folds_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_num_folds_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_dna_level_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_dna_level_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate folding statistics
% this is old, needs to be looked at
function button_fold_stat_Callback(hObject, eventdata, handles)

    folding_count = zeros(size(handles.event_info,1),1);
    level_seq = zeros(size(handles.event_info,1),1);
    DNAdelta = str2num(get(handles.edit_dna_level,'String'));

    round_up_limit = str2num(get(handles.edit_dna_round_up,'String')); % 
    round_down_limit = str2num(get(handles.edit_dna_round_down,'String')); % 

    if isfield(handles,'rawtrace_vect') % only do this if rawtrace exists
        for i=1:size(handles.event_info,1) % save sorted rawtrace
            eventpoints = [(handles.event_info(i,8)) (handles.event_info(i,9))]; % event start stop

            eventtrace_baseline = handles.rawtrace_vect(eventpoints(1):eventpoints(2)) - handles.event_info(i,5);
            fold_count = 0;
            eventlength = length(eventtrace_baseline);

            event_levels_norm_raw = eventtrace_baseline/(-1.*DNAdelta);
            event_levels_norm_floor = floor(event_levels_norm_raw);
            event_levels_norm_rest = event_levels_norm_raw - event_levels_norm_floor;
            event_levels_norm = round(event_levels_norm_raw);

            event_levels_norm_full = zeros(length(event_levels_norm)+1,1);
            event_levels_norm_full(2:end,1) = event_levels_norm;
            num_levels = length(event_levels_norm);



            for j=1:num_levels
                if event_levels_norm_rest(j) > round_up_limit
                    event_levels_norm_full(j+1) = event_levels_norm_floor(j) + 1;
                elseif event_levels_norm_rest(1) < round_down_limit
                    event_levels_norm_full(j+1) = event_levels_norm_floor(j);
                else
                    event_levels_norm_full(j+1) = event_levels_norm_full(j);
                end
                current_diff = event_levels_norm_full(j+1)-event_levels_norm_full(j);
                if current_diff >= 1
                    fold_count = fold_count + current_diff;
                end

            end
            fold_count = fold_count-1;
            folding_count(i) = fold_count;

            %%%%%%%%%%%%
            min_num_of_points = 10;

            % this part counts how many subsequent repetitions there are
            repeating_digit_count = 0;
            current_count_pos = 1;
            for k=1:length(event_levels_norm_full)-1

                if event_levels_norm_full(k) == event_levels_norm_full(k+1)
                    repeating_digit_count(current_count_pos) = repeating_digit_count(current_count_pos) + 1;
                else
                    current_count_pos = current_count_pos +1;
                    repeating_digit_count(current_count_pos) = 0;
                end

            end

            total_digit_count = repeating_digit_count +ones(size(repeating_digit_count));

            % this part removes anything that doesnt repeat at least min_num_of_points
            % implements a minimum level duration
            new_index_pos = 1;
            event_levels_norm_full_mod = 0;
            current_mod_index=1;
            for ii=1:length(total_digit_count)
                if total_digit_count(ii) < min_num_of_points

                else
                    event_levels_norm_full_mod(current_mod_index:current_mod_index+total_digit_count(ii)-1) = ...
                        event_levels_norm_full(new_index_pos:new_index_pos+total_digit_count(ii)-1);
                    current_mod_index = current_mod_index + total_digit_count(ii);

                end
                new_index_pos = new_index_pos + total_digit_count(ii);
            end

            % add zero at end
            event_levels_norm_full_mod(end+1) = 0;

            for k=1:length(event_levels_norm_full_mod)-1
                if event_levels_norm_full_mod(k) ~= event_levels_norm_full_mod(k+1)
                    level_seq(i) = level_seq(i)*10+event_levels_norm_full_mod(k);
                end
            end
    %         if event_levels_norm_full(end) ~= event_levels_norm_full(end-1)
    %             level_seq(i) = level_seq(i)*10+event_levels_norm_full(end);
    %         end



        end
        set(handles.text_status,'String','Finished calculating folding counts')
    else
        set(handles.text_status,'String','Rawtrace file not loaded!!')
    end

    folding_matrix = 0;

    for i=0:max(folding_count)
        folding_matrix(i+1)= length(find(folding_count == i));
    end
    folding_matrix = transpose(folding_matrix);

    unique_lev_seq = unique(level_seq);

    for kk = 1:size(unique_lev_seq,1)
        unique_lev_seq_count(kk) = length(find(level_seq == unique_lev_seq(kk)));
    end
    unique_lev_seq_count = transpose(unique_lev_seq_count);

    lev_seq_all = [unique_lev_seq unique_lev_seq_count];

    [Y,I]=sort(lev_seq_all(:,2));
    lev_seq_all_sorted=lev_seq_all(I,:);

    lev_seq_all_sorted_pct = lev_seq_all_sorted;
    lev_seq_all_sorted_pct(:,2) = lev_seq_all_sorted(:,2)./sum(lev_seq_all_sorted(:,2)).*100;


    set(handles.edit_select_fold_min,'String',num2str(min(folding_count)))
    set(handles.edit_select_fold_max,'String',num2str(max(folding_count)))

    handles.folding_count = folding_count;

    guidata(hObject, handles); % update the handles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_show_levels_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_remove_baseline_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_export_R_Callback(hObject, eventdata, handles)

    % Note: There is an exact copy of this code at the end of fun_plot_R_pore 
    averagecounter = 1;
    alleventcount = 1;

    filenumbers = unique(handles.event_more_info_temp(:,1), 'rows'); % total # of trace files detected

    for filecounter=1:length(filenumbers) % for each trace
        currentfile = filenumbers(filecounter);

        eventindexes = find(handles.event_more_info_temp(:,1) == currentfile); % get Rpore value
        
        if ~isempty(eventindexes) % catch strange case of no events

            fileresistance(filecounter) = handles.event_more_info_temp(eventindexes(1),2); % add the Rpore value
            
            fileresistance_from_baseline(filecounter) = handles.voltage/mean(handles.event_info_temp(eventindexes,5)); % average of all event's baselines in that file
            
            filetime(filecounter) = handles.event_more_info_temp(eventindexes(1),4); % time of experiment for that value

        end
    end


    dir_str = strcat(handles.pathname,'plots', filesep); % in the plots dir
    if(~exist(dir_str, 'dir'))
        mkdir(dir_str);
    end
    file_str = strcat(dir_str,'plot_data_Rpore.dat'); % save dat file with Rpore data
    dataoutput = [filetime', fileresistance', fileresistance_from_baseline'];
    dlmwrite(file_str, dataoutput, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_dna_round_up_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_dna_round_up_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_dna_round_down_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_dna_round_down_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_relative_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine which amplitude data to use
function amplitude = fun_get_amplitude(handles)
    voltage = str2double(get(handles.edit_voltage,'String'))/1000; % conductance case
    switch get(handles.pop_amplitude_type,'Value') % type of Amplitude calculation to use
        case 1 % amplitude_minima_maxima
            if get(handles.check_conductance,'Value')
                if get(handles.check_relative,'Value')
                    amplitude = handles.event_info_temp(:,1)./handles.event_info_temp(:,5); % relative conductance case
                else
                    amplitude = handles.event_info_temp(:,1)./voltage; % conductance case
                end
            else
                if get(handles.check_relative,'Value')
                    amplitude = handles.event_info_temp(:,1)./handles.event_info_temp(:,5); % relative case
                else
                    amplitude = handles.event_info_temp(:,1); % current (nA) case
                end
                
            end        
        case 2 % Max Amplitude
            if get(handles.check_conductance,'Value')
                if get(handles.check_relative,'Value')
                    amplitude = handles.event_info_temp(:,2)./handles.event_info_temp(:,5); % relative conductance case
                else
                    amplitude = handles.event_info_temp(:,2)./voltage; % conductance case
                end
            else
                if get(handles.check_relative,'Value')
                    amplitude = handles.event_info_temp(:,2)./handles.event_info_temp(:,5); % relative case
                else
                    amplitude = handles.event_info_temp(:,2); % current (nA) case
                end
            end
        case 3 % Integral / FWHM
            if get(handles.check_conductance,'Value')
                if get(handles.check_relative,'Value')
                    amplitude = handles.event_more_info_temp(:,9)./handles.event_info_temp(:,5); % relative conductance case
                else
                    amplitude = handles.event_more_info_temp(:,9)./voltage; % conductance case
                end
            else
                if get(handles.check_relative,'Value')
                    amplitude = handles.event_more_info_temp(:,9)./handles.event_info_temp(:,5); % relative case
                else
                    amplitude = handles.event_more_info_temp(:,9); % current (nA) case
                end
            end
        otherwise
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_event_location_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_event_location_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_fexptmn_bins_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_fexptmn_bins_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_lexptmn_bins_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_lexptmn_bins_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_firstextrapointmeanmax_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_firstextrapointmeanmax_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_lastextrapointmeanmax_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_select_lastextrapointmeanmax_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_reset_extrapointmeans_Callback(hObject, eventdata, handles)

    set(handles.text_status,'String','Reseting Extra Point Mean Max Values')
    set(handles.edit_select_firstextrapointmeanmax,'String',num2str(max(handles.events_minmax(:,8))))
    set(handles.edit_select_lastextrapointmeanmax,'String',num2str(max(handles.events_minmax(:,10))))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% radiobutton first extra point mean press
function radio_2_firstextrapointmean_Callback(hObject, eventdata, handles)

    % update other radio button's states
    set(handles.text_status,'String','Plotting First Extra Point Mean Hist')
    set(handles.radio_2_scatter,'Value',0);
    set(handles.radio_2_base,'Value',0);
    set(handles.radio_2_dwell,'Value',0);
    set(handles.radio_2_amplitude,'Value',0);
    set(handles.radio_2_current,'Value',0);
    set(handles.radio_2_R_pore,'Value',0);
    set(handles.radio_2_scatter_integral,'Value',0);
    set(handles.radio_2_integral_hist,'Value',0);
    set(handles.radio_2_scatter_heat_map,'Value',0);
    set(handles.radio_2_firstextrapointmean,'Value',1);
    set(handles.radio_2_lastextrapointmean,'Value',0);

    fun_plot_first_extrapoints_hist(handles,'X'); % call scatter int plot function
    guidata(hObject, handles); % update handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots first extrapoint mean histogram
function fun_plot_first_extrapoints_hist(handles,location)
    
    colormap('default')
    
    if strcmp(location,'X') % plot it here in the GUI
        cla(handles.axes_graphs,'reset')
        axes(handles.axes_graphs)
    end

    if strcmp(location,'n') % plot it in a new figure
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        figure()
    end

    bins = str2num(get(handles.edit_fexptmn_bins,'String')); % number of bins


    [n,xout] = hist(handles.events_minmax_temp(:,8),bins); % plot it
    hist(handles.events_minmax_temp(:,8),bins)
    binwidth = xout(2)-xout(1);
    set(handles.edit_diff_fptmn_bl,'String',num2str(mean(handles.events_minmax_temp(:,8))))
    set(handles.edit_diff_fptmn_bl_STD,'String',num2str(std(handles.events_minmax_temp(:,8))))
    
    set(handles.edit_2_bin_width,'String',num2str(binwidth))

    hXLabel = xlabel('Mean Value (nA)');
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

    if strcmp(location,'n') % set correct limits
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
    end
 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_2_lastextrapointmean_Callback(hObject, eventdata, handles)

    % update other radio button's states
    set(handles.text_status,'String','Plotting Last Extra Point Mean Hist')
    set(handles.radio_2_scatter,'Value',0);
    set(handles.radio_2_base,'Value',0);
    set(handles.radio_2_dwell,'Value',0);
    set(handles.radio_2_amplitude,'Value',0);
    set(handles.radio_2_current,'Value',0);
    set(handles.radio_2_R_pore,'Value',0);
    set(handles.radio_2_scatter_integral,'Value',0);
    set(handles.radio_2_integral_hist,'Value',0);
    set(handles.radio_2_scatter_heat_map,'Value',0);
    set(handles.radio_2_firstextrapointmean,'Value',0);
    set(handles.radio_2_lastextrapointmean,'Value',1);

    fun_plot_last_extrapoints_hist(handles,'X'); % call scatter int plot function
    guidata(hObject, handles); % update handles
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots first extrapoint mean histogram
function fun_plot_last_extrapoints_hist(handles,location) 
   
    colormap('default')
    
    if strcmp(location,'X') % plot it here in the GUI
        cla(handles.axes_graphs,'reset')
        axes(handles.axes_graphs)
    end

    if strcmp(location,'n') % plot it in a new figure
        ylimits = get(handles.axes_graphs,'YLim');
        xlimits = get(handles.axes_graphs,'XLim');
        figure()
    end

    bins = str2num(get(handles.edit_lexptmn_bins,'String')); % number of bins


    [n,xout] = hist(handles.events_minmax_temp(:,10),bins); % plot it
    hist(handles.events_minmax_temp(:,10),bins)
    binwidth = xout(2)-xout(1);
    set(handles.edit_2_bin_width,'String',num2str(binwidth))

    hXLabel = xlabel('Mean Value (nA)');
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

    if strcmp(location,'n') % set correct limits
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
    end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_diff_fptmn_bl_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_diff_fptmn_bl_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_diff_fptmn_bl_STD_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_diff_fptmn_bl_STD_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
