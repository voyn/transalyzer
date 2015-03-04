function varargout = GUI_detect(varargin)
% GUI_DETECT M-file for GUI_detect.fig
%      GUI_DETECT, by itself, creates a new GUI_DETECT or raises the existing
%      singleton*.
%
%      H = GUI_DETECT returns the handle to a new GUI_DETECT or the handle to
%      the existing singleton*.
%
%      GUI_DETECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_DETECT.M with the given input arguments.
%
%      GUI_DETECT('Property','Value',...) creates a new GUI_DETECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_detect_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_detect_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help GUI_detect

% Last Modified by GUIDE v2.5 02-Dec-2013 17:47:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_detect_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_detect_OutputFcn, ...
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
function GUI_detect_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to GUI_detect (see VARARGIN)

    % Choose default command line output for GUI_detect
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    if nargin == 3,
        initial_dir = pwd; % current directory is initial directory
    elseif nargin > 4 % otherwise try to open specified folder
        if strcmpi(varargin{1},'dir')
            if exist(varargin{2},'dir')
                initial_dir = varargin{2};
            else
                errordlg({'Input argument must be a valid',...
                         'folder'},'Input Argument Error!')
                return
            end
        else
            errordlg('Unrecognized input argument',...
                     'Input Argument Error!');
            return;
        end
    end

    handles.currentsegment = 0; % This is for loading ABF files
    handles.total_time_so_far = 0; % experimental time tracking

    % Populate the listbox with files in current dir
    fun_load_listbox(initial_dir,handles)

    % display the logo
    logo=importdata('transalyzer.png');
    axes(handles.axes_trace)
    cla
    image(logo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fun_load_listbox
% Populate the listbox with files in current directory.
% Inputs:
% dir_path - the string of the path to load
%          - when initializing pwd (current Matlab directory) is used
function fun_load_listbox(dir_path, handles) % 

    cd (dir_path) % go to the folder
    dir_struct = dir(dir_path); % list files and dirs

    handles.dir_path = dir_path;

    [sorted_names,sorted_index] = sortrows({dir_struct.name}');
    handles.file_names = sorted_names;
    handles.is_dir = [dir_struct.isdir];
    handles.sorted_index = sorted_index;

    guidata(handles.GUI_detect_fig,handles)

    set(handles.listbox_main,'String',handles.file_names,...
        'Value',1) % set to one otherwise, listbox may dissapear

    set(handles.text_path,'String',pwd) % update path


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = GUI_detect_OutputFcn(hObject, eventdata, handles) 

    % Get default command line output from handles structure
    varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_exit_Callback(hObject, eventdata, handles)
    % close the GUI
    delete(gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%              listbox_main_Callback             %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executes on selection change in listbox_main.
% This is used to browse through directories or open files.
% If file is of unknown type it will assume it is a trace file and pass its
% path to fun_load_trace_file
% the correct file type must be selected by user in the dropdown menu!
function listbox_main_Callback(hObject, eventdata, handles) % when clicking on a file in the listbox

    get(handles.GUI_detect_fig,'SelectionType');
    
    % If double click
    if strcmp(get(handles.GUI_detect_fig,'SelectionType'),'open')
        
        % find what file or folder was clicked
        index_selected = get(handles.listbox_main,'Value');
        
        % get the list of files and folders
        file_list = get(handles.listbox_main,'String');
        
        % Item selected in list box
        filename = file_list{index_selected};
        
        % If folder
        if  isdir(filename)
            
            cd (filename);
            
            % Load list box with new folder.
            fun_load_listbox(pwd,handles);
            handles = guidata(gcbo);
            
        else % if file
            
            % set the file number
            % note first two 'files' ('.', '..') are skipped, since they are the 
            % representations for the current and parent directory 
            set(handles.edit_file_number,'String',num2str(index_selected-2)) 
            
            % always read 1st segment when opening new ABF/TDMS files
            handles.currentsegment = 0; 
            
            % split into path, filename, and extension
            [pathva, name, ext] = fileparts(filename); 
            
            % variable to track if file loading worked
            handles.file_load_successful = 0; 
            
            % cases for different file extensions, add more later
            switch ext 
                
                case '.fig' % Matlab figure
                    
                    % Open FIG-file with guide command.
                    guide (filename)
                    
                    % Use open for other file types.
                    
                otherwise % assume this must be a trace file
                    
                    % try to load it
                    fun_load_trace_file(filename, handles.text_path, handles.currentsegment) 
                    
                    % get updated handles
                    handles=guidata(gcbo);
                    
                    if handles.file_load_successful == 1 % if it was loaded, calculate resistence and plot the trace in GUI
                        
                        fun_calculate_resistance()
                        set(handles.check_highlight_events,'Value',0); % events not detected yet
                        
                        fun_trace_plot('here',0) % plot trace in GUI
                    end
            end
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function listbox_main_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%              Load the File Function            %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function fun_load_trace_file(filename, pathva, currentsegment)
% This function will try to open the trace file and filter it if necessary.
% 
% It will also populate the drop down box for the number of segments.
%
% Inputs:
% filename - complete path to trace file
% currentsegment - which segment of the trace to open
%
% Outputs:
% handles.file_load_successful - file opened correctly 1, or not 0
% handles.trace - current value of the trace in nA
% handles.time_vector - timepoints for each current point
% handles.timestep - time between each point
% handles.total_number_of_segments - total number of segments in this trace
% 
% handles.filter_delay_time - not fully implemented yet, delay introduced
% by filter
% handles.filtered_trace = trace after filtering
% handles.trace_runavg_fwd - moving average
% handles.sigma_window - the STD of the trace
%
function fun_load_trace_file(filename, ~, currentsegment)

                % get handles
                handles=guidata(gcbo); 
                
                % to track if loading worked
                handles.file_load_successful = 0; 
                
                % update handles
                guidata(gcbo, handles);
                
                % get the trace file type
                file_type = get(handles.pop_file_type,'Value');
                
                switch  file_type % switch based on type of trace
                    
                    case 1 % LabView Datalog Binary File
                                                
                        % load the file
                        [handles.trace, handles.time_vector, handles.timestep, handles.file_load_successful] = readlabviewbinaries(filename);
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % this is ugly fix for the odd record size problem
                        % in some LabView binary files
                        handles.trace(handles.trace > 1E20) = [];
                        handles.trace(handles.trace < -1E20) = [];
                        handles.trace(abs(handles.trace) < 1E-20) = [];
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        % generate the time vector
                        handles.time_vector = handles.timestep:handles.timestep:handles.timestep*length(handles.trace);
                        
                        % Apply alpha beta values here
                        alphabeta = str2double(get(handles.edit_ab_gain_value,'String'));
                        handles.trace = handles.trace./alphabeta;
                        
                    case 2 % ATF file (Axon Text File)
                        
                        % try not to use this, it's slow
                        
                        fid = fopen(filename);
                        if fid > 0 % if it was opened
                            for i = 1:10 % lines in header to be skipped
                                line = fgetl(fid);
                            end
                            
                            % read values
                            dataread = fscanf(fid, '%f\t%f', [2 inf]);
                            
                            handles.timestep = dataread(1,3)-dataread(1,2); % calculate timestep
                            handles.trace = dataread(2,:)./1000; % switch to nA
                                                  
                            % be careful at this next part, we are rewriting
                            % experimental time! 
                            % let's say you digitize at 500kSamples/second
                            % you would expect handles.timestep to be 2 us
                            % but it may be represented as
                            % 1.999999999999999999 us in Matlab, take care
                            
                            % create new time vector based on time step,
                            % this closes in gaps, which may be present in
                            % the time vector dataread(1,:);
                            handles.time_vector = handles.timestep:handles.timestep:handles.timestep*length(handles.time_vector);
                            
                            handles.file_load_successful = 1; % worked
                        else
                            handles.file_load_successful = 0; % :(
                        end
                        
                        fclose(fid); % close it up
                        
                    case 3 % Text file with two columns, same as above but without header
                        
                        % try not to use this, it's slow
                        
                        fid = fopen(filename); % open file
                        if fid > 0
                            
                            % read values
                            dataread = fscanf(fid, '%f\t%f', [2 inf]);
                            
                            handles.timestep = dataread(1,3)-dataread(1,2); %% column 1 is time
                            handles.trace = dataread(2,:); %% column 2 is current
                            
                            % handles.time_vector = dataread(1,:); % time vector
                            
                            % be careful at this next part, we are rewriting
                            % experimental time!
                            % let's say you digitize at 500kSamples/second
                            % you would expect handles.timestep to be 2 us
                            % but it may be represented as
                            % 1.999999999999999999 us in Matlab, take care
                            
                            % create new time vector based on time step,
                            % this closes in gaps, which may be present in
                            % the time vector dataread(1,:);
                            handles.time_vector = handles.timestep:handles.timestep:handles.timestep*length(handles.trace);
                            
                            handles.file_load_successful = 1; % worked
                        else
                            handles.file_load_successful = 0; % :(
                        end
                        
                        fclose(fid); % close it up
                        
                    case 4 % no time, current file, read time from the Timestep text box and generate time vector with it
                        
                        % try not to use this, it's slow
                        
                        fid = fopen(filename); % open file
                        
                        if fid > 0
                            dataread = fscanf(fid, '%f', [1 inf]); % read data
                            
                            handles.timestep = str2double(get(handles.edit_timestep,'String'))/1E6; % get timestep from text box
                            handles.trace = dataread(1,:); %% column 1 is current
                            handles.time_vector = handles.timestep:handles.timestep:handles.timestep*length(handles.trace);
                            
                            handles.file_load_successful = 1; % worked
                        else
                            handles.file_load_successful = 0; % :(
                        end
                        
                        fclose(fid); % close it up
                        
                    case 5 % ABF file
                        
                        % binary files (like ABF) are fast
                        
                        segment_size_in_sec = str2double(get(handles.edit_segment_size_in_sec,'String')); % how big each segment is in seconds
                        
                        if currentsegment == 0 % read segment 1 (first 3 sec)
                            
                            [handles.trace,si,lActualAcqLength]=abfload(filename,'start',0,'stop',segment_size_in_sec);
                            set(handles.edit_segment_number,'String','1'); % update segment number textbox
                            
                        else % read some other segment #
                            
                            [handles.trace,si,lActualAcqLength]=abfload(filename,'start',segment_size_in_sec*(currentsegment-1),'stop',segment_size_in_sec*currentsegment);
                            set(handles.edit_segment_number,'String',num2str(currentsegment)); % update segment number textbox
                            
                        end
                        
                        %%% NOTE!!!
                        % sometimes the timestep is not quite exact
                        % instead of 2 us it is 2.000000000001 us
                        % be careful here
                        handles.timestep = si*10^(-6); % convert to sec
                        
                        handles.trace = transpose(handles.trace/1000); % convert to nA, and switch to correct form
                        handles.time_vector = handles.timestep:handles.timestep:handles.timestep*length(handles.trace); % generate time vector
                        
                        numpoints_in_segment = round(segment_size_in_sec/(handles.timestep)); % number of points in a segment
                        handles.total_number_of_segments = int32(idivide(int32(lActualAcqLength),int32(numpoints_in_segment))); % total number of segments in this ABF file
                        
                        for iii=1:handles.total_number_of_segments % generate a cell array to be used in the segment selection popup menu
                          segment_vector{iii} = num2str(iii);
                        end
                        
                        segment_vector = transpose(segment_vector); % flip it
                        set(handles.pop_abf_segment_number,'String',segment_vector) % generate new segment selection popup menu
                        
                        
                        if length(handles.trace) == numpoints_in_segment
                            
                            handles.file_load_successful = 1; % worked
                            
                            if currentsegment == 0 % update current segment number
                                
                                set(handles.pop_abf_segment_number,'Value',1)
                                currentsegment = 1; % current segment is 1
                                
                            else
                                
                                set(handles.pop_abf_segment_number,'Value',currentsegment) % update popup menu
                                
                            end
                        else
                            
                            handles.file_load_successful = 0; % :(
                        end
                        
                    case 6 % TDMS
                        
                        % custom settings may need to be modified for your
                        % particular case
                        segment_size_in_sec = 0.4; % if this is too big it will crash computer
                        Fs = 5e6; % aquire at 5MHz
                        
                        set(handles.edit_segment_size_in_sec,'String',num2str(segment_size_in_sec)) % set the segment size in the text box
                        handles.timestep = 1./Fs; % 
                        numpoints_in_segment = segment_size_in_sec*Fs; % number of points in a segment
                        
                        if currentsegment == 0 % read segment 1 (first 3 sec)
                            
                            converted_tdms_data = TDMS_readTDMSFile(filename,'SUBSET_GET',[1 numpoints_in_segment],'SUBSET_IS_LENGTH',true);
                            set(handles.edit_segment_number,'String','1'); % update segment number textbox
                            
                        else % read segment #    
                            
                            converted_tdms_data = TDMS_readTDMSFile(filename,'SUBSET_GET',[numpoints_in_segment*(currentsegment-1)+1 numpoints_in_segment],'SUBSET_IS_LENGTH',true);
                            set(handles.edit_segment_number,'String',num2str(currentsegment)); % update segment number textbox
                            
                        end
                        
                        [~,metaStruct] = TDMS_readTDMSFile(filename,'GET_DATA_OPTION','getnone'); % get it's properties to see how long it is

                        handles.total_number_of_segments = floor(metaStruct.numberDataPoints(3)/numpoints_in_segment); % total number of segments in this TDMS file

                        for iii=1:handles.total_number_of_segments % generate a cell array to be used in the segment selection popup menu
                          segment_vector{iii} = num2str(iii);
                        end
                        
                        segment_vector = transpose(segment_vector); % flip it
                        set(handles.pop_abf_segment_number,'String',segment_vector) % generate new segment selection popup menu
                        
                        
                        % NOTE: this next part should be changed to suit
                        % your own needs depending on how you recorded your
                        % current signal
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Written by Chin-hsuan (Jennifer) Chen
                        % http://www.mathworks.com/matlabcentral/fileexchange/30023-tdms-reader
                        % include version_2p5_Final
                        % Please change filename, resistor_fb, Fs (Sampling frequency) accordingly
                        % Total Translocation Current: data_tot_current
                        % Baseline Current: baseline_current
                        % Translocation current without Baseline: dna_recover_current_wo_baseline
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        
                        % Coefficient for converting the recorded voltage back to current 
                        const1 = 0;
                        const2 = 2.9883e-24;
                        const3 = 32.9435;

                        Fs = 5e6;
                        resistor_fb = 100e6;

                        % extract the data from TDMS Format
                        voltage_signal = converted_tdms_data.data{:,4};
                        handles.baseline_current = -converted_tdms_data.data{:,3}/resistor_fb;
                        
                        
                        dna_recover_current_wo_baseline = transpose(const1 + const2.*exp(voltage_signal(:).*const3));

                        handles.trace = 1E9*(handles.baseline_current - dna_recover_current_wo_baseline);
                        handles.time_vector = handles.timestep:handles.timestep:handles.timestep*length(handles.baseline_current); % generate time vector
%                         assignin('base', 'voltage_signal', voltage_signal)
%                         assignin('base', 'baseline_current', handles.baseline_current)
%                         assignin('base', 'baseline_current_mean', mean(handles.baseline_current))
%                         assignin('base', 'dna_recover_current_wo_baseline', dna_recover_current_wo_baseline)
%                         assignin('base', 'dna_recover_current_wo_baseline_mean', mean(dna_recover_current_wo_baseline))
%                         assignin('base', 'trace', handles.trace)

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        if length(handles.trace) == numpoints_in_segment
                            
                            handles.file_load_successful = 1; % worked
                            
                            if currentsegment == 0 % update current segment number
                                
                                set(handles.pop_abf_segment_number,'Value',1)
                                currentsegment = 1; % current segment is 1
                                
                            else
                                
                                set(handles.pop_abf_segment_number,'Value',currentsegment) % update popup menu
                                
                            end
                        else
                            handles.file_load_successful = 0; % :(
                        end
                         
                    case 7 % LabView Datalog Binary File (current only)
                        
                        % reading in DTLG files is a custom implementation
                        % depending on LabView script used
                        % the DTLG format is not open, this implementation
                        % was created by analysis of DTLG files generated
                        % by a particular scipt using hexeditors
                        
                        % see
                        % http://galileo.phys.virginia.edu/~pmm3w/misc/file_format.txt
                        
                        % TODO generate a list of all record positions in
                        % the file, useful for long files
%                         if exist(strcat(filename,'_record_positions'),'file') == 0
%                             [~, ~, ~] = readlabviewbinaries_generate_position_file(filename);
%                         end
                        
                        % load the file
                        [handles.trace, handles.time_vector, handles.timestep, handles.file_load_successful] = readlabviewbinaries_readall(filename);
                        
                        if size(handles.trace,1) > size(handles.trace,2)
                            handles.trace = transpose(handles.trace);
                        end      

%                       handles.trace(handles.trace > 1E20) = [];
%                       handles.trace(handles.trace < -1E20) = [];
%                       handles.trace(handles.trace < 1E-20) = [];
                        handles.time_vector = handles.timestep:handles.timestep:handles.timestep*length(handles.trace);
                        
                        % Apply alpha beta values here
                        alphabeta = str2double(get(handles.edit_ab_gain_value,'String'));
                        handles.trace = handles.trace./alphabeta;
                        
                        
                        
                        % NOTE: if you want to add new file formats, they should go here as a new case statement %
                        
                    otherwise
                end
                
                if handles.file_load_successful == 0
                    
                    set(handles.edit_status,'String','File Loading Error')
                    
                else % go here if the trace was loaded successfully
                    
                    % display message say file loaded successfully
                    if get(handles.pop_file_type,'Value') == 5 % if its an ABF give the segment too
                        set(handles.edit_status,'String',strcat(filename, ' loaded. Segment ', {' '},num2str(currentsegment)))
                    else
                        set(handles.edit_status,'String',strcat(filename, ' loaded.'))
                    end
                    
                    % update path string
                    set(handles.text_path,'String',strcat(filename))
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%% Filtering Code begins HERE
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    filter_is_on = get(handles.check_use_filter,'Value');
                    handles.filter_delay_time = 0;
                    
                    if filter_is_on == 1 %% if filtering is turned on
                        
                        % the cutoff frequency
                        handles.cutoff_frequency = str2double(get(handles.edit_filter_fc,'String')); 
                        
                        % code for different types of filters
                        switch get(handles.pop_filter_type,'Value')   
                            case 1 % Gaussian filter
                                
                                % Modified from:
                                
                                % http://www.mathworks.com/matlabcentral/fi
                                % leexchange/12606-1d-gaussian-lowpass-filter
                                    
                                b_gauss=0;
                                a_gauss=3.011*handles.cutoff_frequency;
                                N_gauss=ceil((1/sqrt(2*pi))/(handles.timestep*handles.cutoff_frequency));   
                                %filter half-width, excluding midpoint
                                %Width N corresponds to at least +-3 sigma which captures at least 99.75%
                                %of area under Normal density function. sigma=1/(a*sqrt(2pi)).
                                
                                L_gauss = 2*N_gauss+1;                %full length of FIR filter
                                for k_gauss = -N_gauss:N_gauss
                                    b_gauss(k_gauss+N_gauss+1)=3.011*(handles.cutoff_frequency*handles.timestep)*exp(-pi*(a_gauss*k_gauss*handles.timestep)^2);
                                end;
                                
                                %b(k) coeffs computed above will add to almost exactly unity, but not 
                                %quite exact due to finite sampling and truncation at +- 3 sigma.  
                                %Next line adjusts to make coeffs b(k) sum to exactly unity.
                                filter_b = b_gauss/sum(b_gauss);
                                filter_a = 1;
                                
                                handles.filter_delay_time = N_gauss*handles.timestep; % delay in seconds
                                
                            case 2 % Butterworth 2nd order
                                handles.filter_order = 2;
                                [filter_b,filter_a] = butter(handles.filter_order,handles.cutoff_frequency/(1/(2*handles.timestep)),'low');
                                
                                handles.filter_delay_time = 0; % delay in seconds (TODO)
                                
                            case 3 % Butterworth 4th order
                                handles.filter_order = 4;
                                [filter_b,filter_a] = butter(handles.filter_order,handles.cutoff_frequency/(1/(2*handles.timestep)),'low');
                                
                                handles.filter_delay_time = 0; % delay in seconds (TODO)
                                
                            otherwise
                        end
                        
                        if file_type == 6 % special code for TDMS files
                            % note: do not use filtfilt, introduces distortions
                            
                            % filter the trace
                            filtered_voltage_signal = filter(filter_b,filter_a,voltage_signal);
                            
                            dna_recover_current_wo_baseline_filtered = const1 + const2.*exp(filtered_voltage_signal.*const3);

                            handles.filtered_trace = 1E9*(handles.baseline_current - dna_recover_current_wo_baseline_filtered);
                            
                        else % normal code, for all other file types
                            
                            % note: do not use filtfilt, introduces distortions
                            
                            % filter the trace
                            handles.filtered_trace = filter(filter_b,filter_a,handles.trace);
                        end
                        
                        
                    else %%% for unfiltered traces
                        
                        % trick to bypass need to split filtered and
                        % unfiltered data
                        
                        handles.filtered_trace = handles.trace; % no filtering
                        
                    end %% end filtering if statement
                    
                    trace_mean = mean(handles.filtered_trace); % calc mean current level for the whole trace
                    
                    handles.windowSize = str2double(get(handles.edit_window_size,'String'));
                    
                    % which method to use to calculate the baseline
                    switch get(handles.pop_baseline_method,'Value')
                        case 1 % moving avg baseline, normal
                    
                            set(handles.edit_manual_baseline_value,'String',num2str(trace_mean)) % update the textbox if manual mode not on
                            
                            % next we calculate the moving average
                            
                            % Fast code
                            % see Writing Fast Matlab code, Pascal
                            % Getreuer, February 2009, page 26
                            % http://www.sal.ufl.edu/NewComers/matlab_optim
                            % ization_2.pdf

                            % downside is that it requires windows size #
                            % of points to converge
                            
                            handles.trace_runavg_fwd = cumsum(handles.filtered_trace)/handles.windowSize;
                            handles.trace_runavg_fwd(handles.windowSize+1:end) = handles.trace_runavg_fwd(handles.windowSize+1:end) - handles.trace_runavg_fwd(1:end-handles.windowSize);
                            
                            % small problem that it takes first windowSize
                            % points to converge, fix it here
                            % run it again backwards
                            trace_first_window_flipped = fliplr(handles.filtered_trace(1:4*handles.windowSize));
                            trace_runavg_temp = cumsum(trace_first_window_flipped)/handles.windowSize;
                            trace_runavg_temp(handles.windowSize+1:end) = trace_runavg_temp(handles.windowSize+1:end) - trace_runavg_temp(1:end-handles.windowSize);

                            % flip it
                            trace_first_window_corrected = fliplr(trace_runavg_temp);

                            % in the starting region replace the forward 
                            % with the backward, not elegant, change in
                            % future TODO
                            handles.trace_runavg_fwd(1:handles.windowSize+400) = trace_first_window_corrected(1:handles.windowSize+400);
                             
         
                        case {2, 3} % Manually set initial basline value in nA
                            
                            % this is useful if the total event times
                            % become comparable to the total trace time
                            % such that the average is significantly
                            % distorted
                            handles.trace_runavg_fwd = ones(size(handles.filtered_trace))*str2double(get(handles.edit_manual_baseline_value,'String'));

                        case 4 % Current Histogram Maximum
                            
                            bins_in_current_hist = str2double(get(handles.edit_current_hist_bins,'String'));
                            
                            % make a histogram of all current values
                            [n_temp,xout_temp] = hist(handles.filtered_trace,bins_in_current_hist);
                            
                            % get the position of the peak
                            max_value_baseline_from_hist = xout_temp(n_temp==max(n_temp));
                            
                            % only take one value
                            max_value_baseline_from_hist = max_value_baseline_from_hist(1);
                            
                            % use this value everywhere
                            handles.trace_runavg_fwd = ones(size(handles.filtered_trace))*max_value_baseline_from_hist;
                            
                        case 5 % Baseline from amplifier (TDMS)
                            
                            % if your amplifier calculates the baseline for
                            % you
                            
                            if isfield(handles,'baseline_current') && get(handles.pop_file_type,'Value') == 6 %TDMS
                                
                                handles.trace_runavg_fwd = handles.baseline_current;
                                
                            else
                                
                                set(handles.pop_baseline_method,'Value',1)
                                set(handles.edit_status,'String','Error: No baseline vector available!!')
                            end

                            
                        otherwise
                            
                    end
                    
                                       
                    % number of points in trace
                    set(handles.edit_std_win_size,'String',num2str(length(handles.trace_runavg_fwd)))

                    % this variable is used as the value of sigma
                    handles.sigma_window = 0;
                    
                    % which methods to use for calculating the standard
                    % deviation of the trace
                    switch get(handles.pop_sigma_method,'Value') 
                        case 1 % trace STD
                                                        
                            if filter_is_on == 1 
                                % TODO this is inelegant:
                                % for the first window skip the first 100 datapoints, because of distortions caused by filter
                                % this was fixed by replacing forward
                                % segment with backward? TODO Check it
                                number_of_points_to_skip_at_start_for_std = 100;
                                handles.sigma_window = std(handles.filtered_trace(number_of_points_to_skip_at_start_for_std:end));
                                
                            else
                                handles.sigma_window = std(handles.filtered_trace);
                            end
                            
                        case 2 % Manual (Initial Guess)
                            
                            % read value from GUI text box
                            handles.sigma_window = str2double(get(handles.edit_sigma,'String'));
                            
                            
                        case 3 % Manual (Locked for all iterations)
                            
                            % read value from GUI text box
                            handles.sigma_window = str2double(get(handles.edit_sigma,'String'));
                            
                        case 4 % Moving STD hist peak, normal
                                             
                            % size of window to calc STD
                            % keep it small
                            moving_STD_window_size = str2double(get(handles.edit_sigma_moving_STD_win_size,'String'));
                            
                            sigma_moving = movingstd(handles.filtered_trace(400:end),moving_STD_window_size,'central');
                            
                            % alter number of bins to achieve desired
                            % resolution set in text box
                            desired_STD_resolution = str2double(get(handles.edit_target_res_nA,'String'));
                            
                            number_of_bins_in_current_hist = ceil(abs(max(sigma_moving)-min(sigma_moving))/desired_STD_resolution);
                            
                            % update text box value
                            set(handles.edit_STD_hist_bins,'String',num2str(number_of_bins_in_current_hist))
                            
                            [n_STD,xout_STD] = hist(sigma_moving,number_of_bins_in_current_hist); % make the histogram
                            
                            handles.sigma_window = xout_STD(n_STD==max(n_STD)); % use max point value
                            handles.sigma_window = handles.sigma_window(1); % sometimes multiple bins with the same max height, only take the first
                          
                        otherwise
                            
                    end
                    % show the sigma of the first window
                    set(handles.edit_sigma,'String',num2str(handles.sigma_window))
                    
                    % update values
                    handles.trace_mean = mean(handles.filtered_trace);
                    handles.peakdetectionfactor = str2double(get(handles.edit_peak_detection_factor,'String'));
                    
                end
                set(handles.edit_timestep,'String',num2str(handles.timestep*1E6))
                guidata(gcbo, handles);
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Plot the Trace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function plots the current trace being analyzed. Various features
% can be added. Events can be shown in red with black dots on the start and
% end points. Diagnostic plot mode shows running average as a magenta line, the
% one sigma level above and below the moving average in green, and the
% detection level in red.
%
%
% Inputs:
% location - 'here' to use the handles.axes_trace plot in the GUI, or 'new'
%               to plot it in a new plot window
% show_events_on_trace - 1 also show events found, or 0 don't show events
%               if show events is turned on than handles.check_highlight_events checkbox
%               must be turned on and the vectors handles.event_time_withoutextra and
%               handles.events_local_trace_withoutextra must contain event trace
%               information, while the vectors handles.event_startstop_time and
%               handles.event_startstop_trace must have the start and stop points for all
%               events. These vectors are generated in the function fun_events_detect.
%
% Outputs:
%
%
function fun_trace_plot(location,show_events_on_trace)

    handles=guidata(gcbo); % get handles
    handles.extra_points = str2double(get(handles.edit_extra_points,'String'));
    
    if strcmp(location, 'here') % plot in the big axes at the bottom
        cla(handles.axes_trace,'reset')
        axes(handles.axes_trace)
    end
    
    ylimits = -1; % set for "old" plot
    
    if strcmp(location, 'new') % plot in new figure
        ylimits = get(handles.axes_trace,'YLim'); % save limits
        xlimits = get(handles.axes_trace,'XLim'); % save limits
        figure()
    end
    
    % if there is a lot of data, it is useful to downsample the plot so it
    % does not take forever
    if get(handles.check_downsample,'Value') % are we downsampling? YES
        downsample_number = str2num(get(handles.edit_downsample,'String')); % downsampling factor
        trace_down = downsample(handles.filtered_trace(handles.extra_points:end-handles.extra_points),downsample_number); % downsample trace
        time_down = downsample(handles.time_vector(handles.extra_points:end-handles.extra_points),downsample_number); % downsample time
        plot(time_down, trace_down, '-b') % plot downsampled trace
    else % No downsampling
        % plot trace but skip extra points at start and end
        plot(handles.time_vector(handles.extra_points:end-handles.extra_points), handles.filtered_trace(handles.extra_points:end-handles.extra_points), '-b')
    end
        
    hold on; % we might want to plot more things...
    
    % diagnostic plot, by default always on, shows running avg, std lines,
    % and detection level
    if get(handles.check_diagnostic_plot,'Value') == 1

        if get(handles.check_downsample,'Value') % are we downsampling? YES

            % downsample the moving average
            run_avg_ds = downsample(handles.trace_runavg_fwd(handles.extra_points:end-handles.extra_points),downsample_number);
            
            % plot the downsampled running average in magenta
            plot(time_down, run_avg_ds, '-m');
            hold on;

            % plot std lines above and below in green
            plot(time_down, run_avg_ds-handles.sigma_window, '-g');
            plot(time_down, run_avg_ds+handles.sigma_window, '-g');

            % plot detection level in red
            plot(time_down, run_avg_ds+handles.peakdetectionfactor*handles.sigma_window, '-r');
            plot(time_down, run_avg_ds-handles.peakdetectionfactor*handles.sigma_window, '-r');
        
        else % No downsampling
            
            % plot the running average in magenta
            plot(handles.time_vector(handles.extra_points:end-handles.extra_points), handles.trace_runavg_fwd(handles.extra_points:end-handles.extra_points), '-m');
            hold on;

            % plot std lines above and below in green
            plot(handles.time_vector(handles.extra_points:end-handles.extra_points), handles.trace_runavg_fwd(handles.extra_points:end-handles.extra_points)-handles.sigma_window, '-g');
            plot(handles.time_vector(handles.extra_points:end-handles.extra_points), handles.trace_runavg_fwd(handles.extra_points:end-handles.extra_points)+handles.sigma_window, '-g');

            % plot detection level in red
            plot(handles.time_vector(handles.extra_points:end-handles.extra_points), handles.trace_runavg_fwd(handles.extra_points:end-handles.extra_points)+handles.peakdetectionfactor*handles.sigma_window, '-r');
            plot(handles.time_vector(handles.extra_points:end-handles.extra_points), handles.trace_runavg_fwd(handles.extra_points:end-handles.extra_points)-handles.peakdetectionfactor*handles.sigma_window, '-r');
        
        end
        
    end

    % plot events on trace, if we requested it and if checkbox is checked
    % and if there is at least one event to show
    if ((show_events_on_trace == 1) || (get(handles.check_highlight_events,'Value')) && (str2num(get(handles.text_number_of_events_found,'String')) > 0))
        % plot the current points for every event as red dots
        plot(handles.event_time_withoutextra, handles.events_local_trace_withoutextra, '.r');
        hold on;
        % add black dots on start and stop points for each event
        plot(handles.event_startstop_time, handles.event_startstop_trace, '.k','MarkerSize',5); 
        hold on;
    end
    
    if ylimits == -1 % old plot
        axis([0 max(handles.time_vector) min(handles.filtered_trace(handles.windowSize:end-handles.extra_points)) max(handles.filtered_trace(handles.windowSize:end-handles.extra_points))])
    else % new figure plot
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)])
    end
    
    xlabel('Time (s)')
    ylabel('Current (nA)')
    zoom on; % turn on figure zooming
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% The Event Detection Function %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the core function of this GUI, it will detect events inside a
% trace using a number of approaches.
%
% Called By:
% function fun_cycle_all - automated event detection when cycling through
%                   files/segments/folders. Called multiple times.
% function button_find_events_Callback - when the Find Events button is
%                   pressed.
%
% Inputs:
% plottingon - show the events found on the trace, this is disabled except
%               when finding manually by clicking the Find Events Button
%
% Outputs:
%
function fun_events_detect(plottingon)
    handles=guidata(gcbo); % get handles
    
    %%%%%
    % clear all temp variables here, garbage keeping
    %%%%%
    if isfield(handles,'event_info')
        handles = rmfield (handles,'event_info');
        handles = rmfield (handles,'event_more_info');
        handles = rmfield (handles,'events_local_trace');
        if isfield(handles,'event_unfiltered_trace')
            handles = rmfield (handles,'event_unfiltered_trace');
        end
        handles = rmfield (handles,'event_time');
        handles = rmfield (handles,'events_local_trace_withoutextra');
        handles = rmfield (handles,'event_time_withoutextra');
        handles = rmfield (handles,'event_startstop_trace');
        handles = rmfield (handles,'event_startstop_time');
        handles = rmfield (handles,'event_start_trace_new');
        handles = rmfield (handles,'event_stop_trace_new');
    end
    
    % how many extra points on either side to take, typically 50 or 100
    handles.extra_points = str2double(get(handles.edit_extra_points,'String'));
    
    % event direction 1 - down, 2 - up, 3 - both
    event_direction = get(handles.pop_event_direction,'Value');
    
    % how many iterations to do
    handles.iterations = str2double(get(handles.edit_number_iterations,'String'));
    
    % need to be at least extra_points away from edge to detect
    % this was increased to 2x after some issues, can't recall details
    handles.skip_first = 2*handles.extra_points; % points to skip at start and end of trace
    handles.skip_last = 2*handles.extra_points; % points to skip at start and end of trace
    
    % minimum number of points between sucessive events
    % note that this will merge together these events with smaller separations, although the more
    % likely scenario is that some brief noise spike has caused you to pass
    % the threshold during the course of a single event
    gapsize = str2double(get(handles.edit_skip_points_between_events,'String'));
    
    % are we filtering?
    filter_is_on = get(handles.check_use_filter,'Value');
         
    % dwell time limits here
    handles.min_dwell = str2double(get(handles.edit_minimum_dwell_time,'String')); % min dwell time in ms
    handles.max_dwell = str2double(get(handles.edit_maximum_dwell_time,'String')); % max dwell time in ms
    handles.max_dwell_points = ceil(handles.max_dwell/1000/handles.timestep); % max dwell time in points
    
    % begin iteration loop
    for current_iteration = 1 : handles.iterations + 1
        
        % iteration loop variables initialization
        handles.filtered_trace_events_replaced_with_baseline = handles.filtered_trace;
        handles.filtered_trace_events_replaced_with_NaN = handles.filtered_trace;
        handles.filtered_trace_events_removed = handles.filtered_trace;
        
        event_counter = 0; % count your events

        % shift baseline to zero
        trace_events_above_zero = handles.filtered_trace - handles.trace_runavg_fwd;

        % any point outside threshold corresponds to an event
        if event_direction == 1 % down events
            event_indexes = find(trace_events_above_zero < -1*handles.peakdetectionfactor*handles.sigma_window);
        elseif event_direction == 2 % up events
            event_indexes = find(trace_events_above_zero > handles.peakdetectionfactor*handles.sigma_window);
        elseif event_direction == 3 % up and down events
            event_indexes = find(abs(trace_events_above_zero) > handles.peakdetectionfactor*handles.sigma_window);
        end

        % now we will use a trick to speed things up, for most events there
        % will be lots of consecutive points corresponding to the time the
        % current level is above the detection level. So subsequent points
        % are from the same event but gaps probably represent the
        % boundaries of each event.
        % We find these gaps by shifting the vector and subtracting it from
        % itself.
        % These gaps then become our initial guesses to find where the
        % baseline was crossed. This approach save time in having FOR loops
        % scan over large regions.
        
        % shift vector with points one over    
        event_indexes_shifted = transpose(circshift(event_indexes',-1));

        % subtract shifted from not shifted
        edge_detect = event_indexes_shifted - event_indexes;

        % If difference is above the minimum time between events, treat
        % them as different events. If they happend to be the same event, this will
        % be caught by some conditional statements later on.
        edge_index_index = find(abs(edge_detect) > abs(1 + gapsize));
        
        % calculate the starting and stopping points where threshold is
        % crossed
        starting_edge_index = event_indexes_shifted(edge_index_index);
        stopping_edge_index = event_indexes(edge_index_index);

        % shift back
        if length(starting_edge_index) > 1 % if more than 1 event, shift start positions to correct place
            starting_edge_index = transpose(circshift(starting_edge_index',1));
        end

        % this is added to make sure you dont detect the same event twice
        % the start of a new event must be after the end of the previous
        % event
        previous_end_point = 1;

        % here we do another trick to find all baseline crossing points
        
        % first shift the trace so that the baseline is zero
        filtered_trace_shifted_fwd = handles.filtered_trace - handles.trace_runavg_fwd;
        
        % convert positive values to 1 negative to -1 and zero values to 0
        filtered_trace_shifted_sign_change_fwd = sign(filtered_trace_shifted_fwd); % if positive = 1 if negative = -1 if zero = 0

        % this is an error trap for the case where we actually have a value
        % exactly equal to zero, TODO not elegant, better way?
        % for now treat zero as positive number
        filtered_trace_shifted_sign_change_fwd(filtered_trace_shifted_sign_change_fwd==0) = 1;

        % get baseline crossover points, by looking for places where value
        % of adjacent numbers changes by 2
        filtered_trace_crossing_points_fwd = find(abs(diff(filtered_trace_shifted_sign_change_fwd))==2);

        % events not allowed too close to the end, this is the limit
        max_index_for_event_end = length(handles.filtered_trace) - handles.skip_last;
        
        for i = 1:length(starting_edge_index) % for all possible events

            end_of_file_flag = 0; % reached the end of the file yet?
            
            % Mother of all IF statements
            % checks that:
            %  1) event starts after the edge of the first window (TODO)
            %  2) events do not overlap
            %  3) events not too close to the end or start
            %  4) there is enough space on either side to save the extrapoints
            if (starting_edge_index(i) > handles.skip_first) &&...
                    (starting_edge_index(i) > previous_end_point) &&...
                    (stopping_edge_index(i) < (length(handles.trace) - handles.skip_last))

                % first guess about where event starts is location after
                % baseline crossing but before crossing detection level
                % find the last baseline crossing point before the
                % detection level was crossed
                event_start = filtered_trace_crossing_points_fwd(find(filtered_trace_crossing_points_fwd < starting_edge_index(i), 1, 'last'));

                if isempty(event_start) % catch errors here
                    end_of_file_flag = 1; % this will cause this 'event' to be discarded
                    event_start = starting_edge_index(i);
                end
                
                if event_start < previous_end_point % make sure no overlapping events
                    end_of_file_flag = 1; % this will cause this 'event' to be discarded
                end
                
                
                % the official local baseline for this event
                baseline = handles.trace_runavg_fwd(event_start);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % this next part is fast but, prone to error because the moving
                % average often changes due to the event itself
                % this would be the code:
                % event_end =  filtered_trace_crossing_points_fwd(find(filtered_trace_crossing_points_fwd > stopping_edge_index(i), 1, 'first'))+1;
                % instead we use this:
                % TODO, this part is ugly, should be recoded with vector operations if possible
                % first guess about where event ends is location before
                % crossing detection level, we will find the exact point
                % later in the code, this is just a guess
                event_end =  stopping_edge_index(i);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % find the end point based on the type of event
                if event_direction == 1 % down event
                    
                    increment_counter = 0;
                    
                    % keep going until the local baseline is crossed
                    while (handles.filtered_trace(event_end) < baseline) && (end_of_file_flag == 0) % has not crossed baseline yet
                    
                        % look for end point, increment each loop
                        event_end=event_end+1;
                        increment_counter = increment_counter + 1;

                        if event_end > max_index_for_event_end || increment_counter > handles.max_dwell_points % dont get too close to the end of the file
                            end_of_file_flag = 1;
                        end
                    end
                    
                elseif event_direction == 2 % up event
                    
                    increment_counter = 0;
                    
                    % keep going until the local baseline is crossed
                    while (handles.filtered_trace(event_end) > baseline) && (end_of_file_flag == 0) % has not crossed baseline yet
                    
                        % look for end point, increment each loop
                        event_end=event_end+1;
                        increment_counter = increment_counter + 1;

                        if event_end > max_index_for_event_end || increment_counter > handles.max_dwell_points % dont get too close to the end of the file
                            end_of_file_flag = 1;
                            % this event is discarded by a subsequent if
                            % statement
                        end
                    end
                    
                elseif event_direction == 3 %both
                    
                    % look for the next baseline crossing, can be in either
                    % direction
                    while (abs(sign(handles.filtered_trace(event_end)- baseline)-sign(handles.filtered_trace(event_end-1)- baseline))~=2) % has not crossed baseline yet

                        % look for end point, increment each loop
                        event_end = event_end + 1;

                        if event_end > max_index_for_event_end % dont get too close to the end of the file
                            end_of_file_flag = 1;
                            break % ugly, but fast
                            % this event is discarded by a subsequent if
                            % statement
                        end
                        
                    end
                    
                    % TODO here we insert some magic
                    % what if we have events which contain both up and down
                    % portions? See DOI: 10.1021/nl301719a
                    % make sure the point after crossing the baseline is greater than the
                    % detection threshold on the other side
                    % Note: this only handles one changeover!!
                    if abs(filtered_trace_shifted_fwd(event_end)) > handles.peakdetectionfactor*handles.sigma_window
                        
                        % if it is a strage event use the next crossing point
                        event_end = filtered_trace_crossing_points_fwd(find(filtered_trace_crossing_points_fwd > event_end, 1, 'first'))+1;
                        display('UP & DOWN!!!')

                    end
                end

                if isempty(event_end) % catch events crossing end of trace
                    end_of_file_flag = 1;
                    event_end = length(handles.trace); 
                end
              
                % using baseline crossing points very innacurate way to calculate dwell time, because of long tails, filtering effects, always
                % use FWHM instead
                event_dwell = (event_end - event_start)*handles.timestep*1000;

                
                % if dewll time is between the min and max allowed values, and
                % the end of the file was not reached, and there is enough
                % space to grab the extra points at the edges
                if (event_dwell > handles.min_dwell) && (event_dwell < handles.max_dwell) && (end_of_file_flag == 0) ...
                        && (event_start > handles.extra_points) && (length(handles.trace_runavg_fwd)-event_end > handles.extra_points)
                    
                    if current_iteration == handles.iterations + 1 % if this is the last iteration analyze the event

                        event_counter = event_counter + 1; % increment event counter

                        event_trace = handles.filtered_trace(event_start:event_end); % get the trace

                        event_time_vect = 0:handles.timestep*1000:(event_end - event_start)*handles.timestep*1000; % build a time vector (ms)

                        event_average = baseline - mean(event_trace); % event average

                        % is this an up or down event
                        up_or_down = (-1*sign(event_average))/2+0.5;

                        % shifted trace
                        removebaseline_trace = event_trace - baseline.*ones(size(event_trace));

                        if up_or_down == 1 % up events

                            %get detection level at crossover point
                            detectionlevel = (handles.trace_runavg_fwd(starting_edge_index(i)) + handles.peakdetectionfactor*handles.sigma_window);

                            [lmval,minindex] = lmax(event_trace, 0); % calc local maxima
                            event_max = max(event_trace) - baseline; % event max

                        else % down events

                            %get detection level at crossover point
                            detectionlevel = (handles.trace_runavg_fwd(starting_edge_index(i)) - handles.peakdetectionfactor*handles.sigma_window);

                            [lmval,minindex] = lmin(event_trace, 0);  % calc local minima
                            event_max = baseline - min(event_trace); % event max

                        end

                        % calc the integral
                        event_int = -1.*trapz(event_time_vect,removebaseline_trace); % in nA*ms

                        if length(minindex) > 1 % more than 1 local minima/maxima point

                            if up_or_down == 1 % up event
                                amplitude_minima_maxima = mean(event_trace(minindex(1):max(minindex))) - baseline; % amplitude calculated only between first and last maxima
                                minmax_index = find(event_trace == max(event_trace));
                                event_trace_norm = (event_trace - baseline) / amplitude_minima_maxima; % for FWHM calc
                            else % down event
                                amplitude_minima_maxima = baseline - mean(event_trace(minindex(1):max(minindex))); % amplitude calculated only between first and last minima
                                minmax_index = find(event_trace == min(event_trace));
                                event_trace_norm = -(event_trace - baseline) / amplitude_minima_maxima; % for FWHM calc
                            end

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % this has problems, tends to give dwell times shorter
                            % than they should be since the last local minima is
                            % not always at the end of the plateau but may be a bit
                            % ahead, use FWHM instead
                            % FWHM_dwell =
                            % (handles.timestep*1000*(minindex(max(length(minindex)))-1));
                            % % dwell time calc similar to Pedone et al.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            minima_length = length(minindex);
                            minima_min = min(minindex);
                            minima_max = max(minindex);

                            % FWHM code
                            % http://www.mathworks.nl/matlabcentral/fileexchange/10590-fwhm

                            % detect starting FWHM point
                            N = length(event_trace_norm);
                            icount = 2; % start at second point and go forward
                            while sign(event_trace_norm(icount)-0.5) == sign(event_trace_norm(icount-1)-0.5) % go forward until you find the crossing
                                icount = icount + 1;
                                if icount == N + 1 %% if you cant find changeover
                                    icount = N; % keep calm, but FreakOut!!
                                    disp('problem with finding starting FWHM point of event #',{' '},num2str(event_counter))
                                    break % and break the loop
                                end
                            end

                            % interpolate the time of the starting FWHM point
                            interp = (0.5-event_trace_norm(icount-1)) / (event_trace_norm(icount)-event_trace_norm(icount-1));
                            tlead = event_time_vect(icount-1) + interp*(event_time_vect(icount)-event_time_vect(icount-1));


                            % must have a + 1 here in case the min/max is on the
                            % first point after the falling edge, this causes dwell
                            % = 0 errors
                            %icountend = minmax_index(length(minmax_index))+1;   %start search for next crossing at last local minima
                            icountend = N;   %start search for next crossing at the end

                            % work your way backwards
                            while ((icountend >= icount) && (sign(event_trace_norm(icountend)-0.5) == sign(event_trace_norm(icountend-1)-0.5)))
                                icountend = icountend - 1;
                                if icountend == 0 %% if you cant find changeover
                                    icountend = 1; % keep calm, but FreakOut!!
                                    disp('problem with finding ending FWHM point of event #',{' '},num2str(event_counter))
                                    break % and break the loop
                                end
                            end
                            if icountend > icount % sanity check
                                Ptype = 1;
                                interp = (0.5-event_trace_norm(icountend-1)) / (event_trace_norm(icountend)-event_trace_norm(icountend-1));
                                ttrail = event_time_vect(icountend-1) + interp*(event_time_vect(icountend)-event_time_vect(icountend-1));

                                FWHM_dwell = ttrail-tlead;% Full Width at Half Maximum in ms

                                % calc the integral based on the points before and
                                % after the FWHM zone
                                % this prevents long tails from messing up the
                                % integral
                                event_int = -1.*trapz(event_time_vect(icount-1:icountend),removebaseline_trace(icount-1:icountend)); % in nA*ms

                            else % big problem caught here
                                Ptype = 2;
                                disp('Step-Like Pulse, no second edge, if you see this, be afraid!! Very afraid!!')

                                FWHM_dwell = event_dwell; % Bad case

                            end

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        else % a single point (or no points), this is a spike

                            if up_or_down == 1 % up event
                                amplitude_minima_maxima = max(event_trace) - baseline; % take the max as the amplitude
                                minmax_index = find(event_trace == max(event_trace));
                            else % down event
                                amplitude_minima_maxima = baseline - min(event_trace); % take the min as the amplitude
                                minmax_index = find(event_trace == min(event_trace));
                            end

                            minima_length = 0;
                            minima_min = 0;
                            minima_max = 0;

                            event_trace_norm = abs(event_trace - baseline) / amplitude_minima_maxima;

                            % detect starting FWHM point
                            N = length(event_trace_norm);
                            icount = 2; % start at second point and go forward
                            while sign(event_trace_norm(icount)-0.5) == sign(event_trace_norm(icount-1)-0.5) % go forward until you find the crossing
                                icount = icount + 1;
                                if icount == N + 1 %% if you cant find changeover
                                    icount = N; % keep calm, but FreakOut!!
                                    disp(strcat('problem with finding starting FWHM point of event #',{' '},num2str(event_counter)))
                                    break % and break the loop
                                end
                            end

                            % interpolate the time of the starting FWHM point
                            interp = (0.5-event_trace_norm(icount-1)) / (event_trace_norm(icount)-event_trace_norm(icount-1));
                            tlead = event_time_vect(icount-1) + interp*(event_time_vect(icount)-event_time_vect(icount-1));

                            % must have a + 1 here in case the min/max is on the
                            % first point after the falling edge, this causes dwell
                            % = 0 errors
                            %icountend = minmax_index(length(minmax_index))+1;   %start search for next crossing at last local minima
                            icountend = N;   %start search for next crossing at the end

                            % work your way backwards
                            while ((icountend >= icount) && (sign(event_trace_norm(icountend)-0.5) == sign(event_trace_norm(icountend-1)-0.5)))
                                icountend = icountend - 1;
                                %assignin('base', 'icount', i)
                                if icountend == 0 %% if you cant find changeover
                                    icountend = 1; % keep calm, but FreakOut!!
                                    disp('problem with finding ending FWHM point of event #',{' '},num2str(event_counter))
                                    break % and break the loop
                                end
                            end
                            if icountend > icount % sanity check
                                Ptype = 1;
                                %disp('In the land of signal spikes.')
                                interp = (0.5-event_trace_norm(icountend-1)) / (event_trace_norm(icountend)-event_trace_norm(icountend-1));
                                ttrail = event_time_vect(icountend-1) + interp*(event_time_vect(icountend)-event_time_vect(icountend-1));

                                FWHM_dwell = ttrail-tlead;% Full Width at Half Maximum in ms

                                %event_trace_FWHM = event_trace(icount-1:icountend);
                                % calc the integral based on the points before and
                                % after the FWHM zone
                                event_int = -1.*trapz(event_time_vect(icount-1:icountend),removebaseline_trace(icount-1:icountend)); % in nA*ms

                            else % big problem
                                Ptype = 2;
                                disp('Step-Like Pulse, no second edge, if you see this, be afraid!! Very afraid!!')

                                FWHM_dwell = event_dwell; % Bad case
                            end

                        end

                        handles.event_info(event_counter,:) = [amplitude_minima_maxima event_max event_int FWHM_dwell baseline detectionlevel up_or_down event_start event_end];
        %               this vector contains:
        %                 1 amplitude_minima_maxima	
        %                 2 maximum	
        %                 3 integral	
        %                 4 FWHM_dwell	
        %                 5 baseline	
        %                 6 detectionlevel
        %                 7 eventtypeUpDown	down =0 up=1
        %                 8 startpoint	
        %                 9 endpoint

                        % the integral method for determining the current
                        % blockade
                        event_amplitude_int_over_FWHM = event_int/FWHM_dwell; 

                        handles.event_more_info(event_counter,:) = [minima_length minima_min minima_max event_dwell event_amplitude_int_over_FWHM];
        %               this vector contains:
                        % 1 number of local min
                        % 2 index first local min
                        % 3 index last local min
                        % 4 dwell time old (between start and stop)
                        % 5 amplitude good (int_over_FWHM)

                        
                        
                    else % this is reached if it is not the last iteration
                         
                        event_counter = event_counter + 1; % increment event counter
                        handles.event_info(event_counter,:) = [0 0 0 0 0 0 0 event_start event_end];
                        
                        % make new trace file for moving average
                        handles.filtered_trace_events_replaced_with_baseline(event_start:event_end) = handles.trace_runavg_fwd(event_start);
                        
                        % make new trace file for sigma
                        handles.filtered_trace_events_replaced_with_NaN(event_start:event_end) = NaN;
                        
                        
                    end
                    
                    previous_end_point = event_end; % save end, for the event overlap check
                    
                end
            end
        end % all possible events analyzed
        
        % this runs at the end of each iteration
        if current_iteration < handles.iterations + 1 % if this is not the last iteration
            

            if event_counter > 0 % if events were found
                
                for i=event_counter:-1:1 % count down through all found events
                    
                    % remove each event from trace                    
                    handles.filtered_trace_events_removed(handles.event_info(i,8):handles.event_info(i,9)) = [];
                    
                end
                
            end

            
            switch get(handles.pop_baseline_method,'Value')
                case {1, 2, 4} % moving avg or Manual Initial or Current Hist Max

                    % recalculate the moving average using the
                    % removed events baseline
                    handles.trace_runavg_fwd = cumsum(handles.filtered_trace_events_replaced_with_baseline)/handles.windowSize;

                    handles.trace_runavg_fwd(handles.windowSize+1:end) = handles.trace_runavg_fwd(handles.windowSize+1:end) - handles.trace_runavg_fwd(1:end-handles.windowSize);

                    % small problem that it takes first windowSize
                    % points to converge, fix it here
                    trace_first_window_flipped = fliplr(handles.filtered_trace_events_replaced_with_baseline(1:4*handles.windowSize));
                    trace_runavg_temp = cumsum(trace_first_window_flipped)/handles.windowSize;
                    trace_runavg_temp(handles.windowSize+1:end) = trace_runavg_temp(handles.windowSize+1:end) - trace_runavg_temp(1:end-handles.windowSize);

                    trace_first_window_corrected = fliplr(trace_runavg_temp);

                    handles.trace_runavg_fwd(1:handles.windowSize+400) = trace_first_window_corrected(1:handles.windowSize+400);

                case 3 % Manual Locked

                    % nothing

                otherwise
            end

            switch get(handles.pop_sigma_method,'Value')
                case 1 % trace STD

                    if filter_is_on == 1 
                        % TODO this is inelegant:
                        % for the first window skip the first 100 datapoints, because of distortions caused by filter
                        number_of_points_to_skip_at_start_for_std = 100;
                        handles.sigma_window = nanstd(handles.filtered_trace_events_replaced_with_NaN(number_of_points_to_skip_at_start_for_std:end));
                    else
                        handles.sigma_window = nanstd(handles.filtered_trace_events_replaced_with_NaN);
                    end
                case 3
                    % nothing

                case {2}
                    if filter_is_on == 1 
                        % TODO this is inelegant:
                        % for the first window skip the first 100 datapoints, because of distortions caused by filter
                        number_of_points_to_skip_at_start_for_std = 100;
                        sigma_window_temp = nanstd(handles.filtered_trace_events_replaced_with_NaN(number_of_points_to_skip_at_start_for_std:end));
                    else
                        sigma_window_temp = nanstd(handles.filtered_trace_events_replaced_with_NaN);
                    end
                    
                    % update the sigma value only if it's converging
                    if sigma_window_temp <= handles.sigma_window
                        handles.sigma_window = sigma_window_temp;
                    end
                    
                case {4} % Moving STD hist peak, normal
                   
                    % recalculate STD each iteration but only update if the
                    % value decreased
                    
                    moving_STD_window_size = str2double(get(handles.edit_sigma_moving_STD_win_size,'String')); % keep it small
                            
                    sigma_moving = movingstd(handles.filtered_trace_events_removed(400:end),moving_STD_window_size,'central');

                    % alter number of bins to achieve desired
                    % resolution
                    desired_STD_resolution = str2double(get(handles.edit_target_res_nA,'String'));

                    number_of_bins_in_current_hist = ceil(abs(max(sigma_moving)-min(sigma_moving))/desired_STD_resolution);

                    set(handles.edit_STD_hist_bins,'String',num2str(number_of_bins_in_current_hist))

                    [n_STD,xout_STD] = hist(sigma_moving,number_of_bins_in_current_hist); % make the histogram

                    sigma_window_temp = xout_STD(n_STD==max(n_STD)); % use max point value
                    sigma_window_temp = sigma_window_temp(1); % sometimes multiple bins with the same max height, only take the first
                    
                    % update the sigma value only if it's converging
                    if sigma_window_temp < handles.sigma_window
                        handles.sigma_window = sigma_window_temp;
                    end
                    
                otherwise
            end
        end
    end
    
    % initialize trace pointers
    eventpointssum = 1;
    eventpointssum_withoutextra = 1;
    
    if event_counter > 0 % if events were found
        for index=1:size(handles.event_info,1) % for every event found
            
            % start stop points with extra points
            eventpoints = handles.event_info(index,9) - handles.event_info(index,8)+1 + 2*handles.extra_points;
            % start stop points without extra points
            eventpoints_withoutextra = handles.event_info(index,9) - handles.event_info(index,8)+1;
            
            % this is the local trace and time of all events concatenated
            handles.events_local_trace(eventpointssum:(eventpointssum + eventpoints-1)) = handles.filtered_trace(handles.event_info(index,8)-handles.extra_points:handles.event_info(index,9)+handles.extra_points);
            handles.event_time(eventpointssum:(eventpointssum + eventpoints-1)) = handles.time_vector(handles.event_info(index,8)-handles.extra_points:handles.event_info(index,9)+handles.extra_points);
            
            % same but without extra points
            handles.events_local_trace_withoutextra(eventpointssum_withoutextra:(eventpointssum_withoutextra + eventpoints_withoutextra-1)) = handles.filtered_trace(handles.event_info(index,8):handles.event_info(index,9));
            handles.event_time_withoutextra(eventpointssum_withoutextra:(eventpointssum_withoutextra + eventpoints_withoutextra-1)) = handles.time_vector(handles.event_info(index,8):handles.event_info(index,9));
            
            % calc the starting and stoping points
            handles.event_start_trace_new(index) = (eventpointssum+handles.extra_points);
            handles.event_stop_trace_new(index) = (eventpointssum+eventpoints-handles.extra_points-1);
            
            if get(handles.check_save_unfiltered_trace,'Value') % save unfiltered trace
                % (TODO) need more code here to take into account the delay
                % introduced by filtering
                delay_in_points = handles.filter_delay_time/handles.timestep;
                handles.event_unfiltered_trace(eventpointssum:(eventpointssum + eventpoints-1)) = handles.trace(handles.event_info(index,8)-handles.extra_points:handles.event_info(index,9)+handles.extra_points);
                % TODO, this should take into account the filter delay, but
                % it also needs to catch cases where the filter delay
                % causes the index to be less than 1
                %               handles.event_unfiltered_trace(eventpointssum:(eventpointssum + eventpoints-1)) = handles.trace(handles.event_info(index,8)-handles.extra_points-delay_in_points:handles.event_info(index,9)+handles.extra_points-delay_in_points);
                
            end
            
            % update the trace pointers
            eventpointssum = eventpointssum + eventpoints;
            eventpointssum_withoutextra = eventpointssum_withoutextra + eventpoints_withoutextra;
            
            % these are used for overlaying the events found on top of the
            % plot of the original trace
            handles.event_startstop_trace((2*index-1):2*index) = [handles.filtered_trace(handles.event_info(index,8)) handles.filtered_trace(handles.event_info(index,9))]; 
            handles.event_startstop_time((2*index-1):2*index) = [handles.time_vector(handles.event_info(index,8)) handles.time_vector(handles.event_info(index,9))];
        end
        

    end
    
    if isfield(handles,'event_info')
        [event_counter cc] = size(handles.event_info); % count the rows
        assignin('base', 'event_info', handles.event_info)
        assignin('base', 'starting_edge_index', starting_edge_index)
        set(handles.text_number_of_events_found,'String',num2str(event_counter)) % update event counter
    else
        set(handles.text_number_of_events_found,'String','0') % :'(
        set(handles.edit_status,'String','No Events Found!')
    end
        
    guidata(gcbo, handles); % update handles
    
    if (event_counter > 0) && plottingon % if plotting is turned on
        fun_trace_plot('here',1) % show the events on the trace
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_window_size_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_window_size_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_peak_detection_factor_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_peak_detection_factor_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_minimum_dwell_time_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_minimum_dwell_time_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_maximum_dwell_time_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_maximum_dwell_time_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_skip_points_between_events_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_skip_points_between_events_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_extra_points_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_extra_points_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_use_filter_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_filter_fc_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_filter_fc_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_ab_gain_value_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_ab_gain_value_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_save_Callback(hObject, eventdata, handles)
    % save it all
    fun_save_detection_results(handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%               fun_save_detection_results                          %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function creates and saves all of the results
function fun_save_detection_results(handles) % save stuff
if isfield(handles,'event_trace_master') % there are events in the master file
    
    % name saving directory
    if get(handles.check_use_filter,'Value') == 1 
        dir_str = strcat(pwd, filesep, '~transalyzed-',num2str(handles.cutoff_frequency),'Hz', filesep);
    else
        dir_str = strcat(pwd, filesep, '~transalyzed-unfiltered', filesep);
    end
    
    if(~exist(dir_str, 'dir')) % make directory if it doesn't exist
        mkdir(dir_str);
        
    end
    
    if get(handles.check_make_OpenNanopore_files,'Value') == 1
        OpenNanopore_dir = strcat(dir_str,'OpenNanopore');
        if(~exist(OpenNanopore_dir, 'dir')) % make directory
            mkdir(OpenNanopore_dir);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fun_save_par_file(handles, dir_str)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % write analysis_parameters.txt file
    file_str = strcat(dir_str,'analysis_parameters.txt');
    fid = fopen(file_str, 'w');
    fprintf(fid, '%s\r\n', '%% Transalyzer Analysis');
    fprintf(fid, '%s\r\n', strcat('%% ', get(handles.text_version,'String')));
    
    fprintf(fid, '%s\r\n', '%% Analysis Date');
    fprintf(fid, '%s\r\n', date);
    
    fprintf(fid, '%s\r\n', '%% Voltage (mV)');
    fprintf(fid, '%s\r\n', get(handles.edit_voltage,'String'));
    
    fprintf(fid, '%s\r\n', '%% Sampling Frequency (Hz)');
    sampl = 1/handles.timestep;
    fprintf(fid, '%s\r\n', num2str(sampl));
    
    fprintf(fid, '%s\r\n', '%% Local Trace Time (s)');
    ttime = handles.timestep*length(handles.filtered_trace);
    fprintf(fid, '%s\r\n', num2str(ttime));
        
    fprintf(fid, '%s\r\n', '%% Peak Det Factor');
    fprintf(fid, '%s\r\n', get(handles.edit_peak_detection_factor,'String'));
    
    fprintf(fid, '%s\r\n', '%% Sigma (nA)');
    fprintf(fid, '%s\r\n', num2str(handles.sigma_window));
    
    fprintf(fid, '%s\r\n', '%% Extra Points (points)');
    fprintf(fid, '%s\r\n', get(handles.edit_extra_points,'String'));
    
    fprintf(fid, '%s\r\n', '%% Total Events Found');
    fprintf(fid, '%s\r\n', get(handles.text_total_number_of_events_found,'String'));
        
    fprintf(fid, '%s\r\n', '%% Total Time (sec)');
    fprintf(fid, '%s\r\n', num2str(handles.total_time_so_far));
    
    fclose(fid); % close file
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % write analysis_events file
    file_str = strcat(dir_str,'analysis_events');
    dlmwrite(file_str, handles.event_info_master, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % convert it into strings
    extraincell = cellstr(num2str(handles.event_more_info_master));
    
    % concatenate arrays together
    allextrastuffconcatenated = horzcat(handles.event_master_trace_info(:,2:5),extraincell);
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
    
    
    % write analysis_extra_info.txt file
    file_str = strcat(dir_str,'analysis_extra_info.txt');
    fid = fopen(file_str, 'w'); % open file for writing
    [M,N] = size(allextrastuffconcatenated); % get array size
	for i=1:M % this is slow, but works
        for j=1:N
            fprintf(fid, '%s\t',allextrastuffconcatenated{i,j}); % write data point
        end
        fprintf(fid, '\r\n'); % new line
	end
    fclose(fid); % close it
    
    % write analysis_files_info.txt file
    file_str = strcat(dir_str,'analysis_files_info.txt');
    fid = fopen(file_str, 'w'); % open file for writing
    filenamesvector = handles.event_master_trace_info(:,1:2);
    % contains names of files and file number
    %filenamesvector = uniqueRowsCA(filenamesvector, 'rows');
%    filenamesvector = sortrows(filenamesvectortemp,2);
    
    [M,N] = size(filenamesvector); % get array size
	for i=1:M % this is slow, but works
        for j=1:N
            fprintf(fid, '%s\t',filenamesvector{i,j}); % write data point
        end
        fprintf(fid, '\r\n'); % new line
	end
    fclose(fid); % close it
    
    % write analysis_rawtrace file
    file_str = strcat(dir_str,'analysis_rawtrace');
    dlmwrite(file_str, handles.event_trace_master, 'delimiter', '\n','precision', '%.8f','newline', 'pc');
    
    if get(handles.check_save_unfiltered_trace,'Value') % save unfiltered trace
        % write analysis_unfiltered_rawtrace file
        file_str = strcat(dir_str,'analysis_unfiltered_rawtrace');
        dlmwrite(file_str, handles.event_unfiltered_trace_master, 'delimiter', '\n','precision', '%.8f','newline', 'pc');
    end
    
    if get(handles.check_make_OpenNanopore_files,'Value') == 1
        fun_make_new_trace(handles, dir_str, 1) % use 1 for Open Nanopore
        % write OpenNanopore_parameters.txt file
        file_str = strcat(dir_str,'OpenNanopore', filesep,'OpenNanopore_parameters.txt');
        fid = fopen(file_str, 'w');
        fprintf(fid, '%s\r\n', '%% OpenNanopore Parameters');
        fprintf(fid, '%s\r\n', '%% Delta');
        fprintf(fid, '%s\r\n', '? nA');
        fprintf(fid, '%s\r\n', '%% Sigma (nA)');
        fprintf(fid, '%s\r\n', num2str(handles.sigma_window));
        fprintf(fid, '%s\r\n', '%% hBook');
        fprintf(fid, '%s\r\n', '?');
        fprintf(fid, '%s\r\n', '%% Sampling Frequency (Hz)');
        fprintf(fid, '%s\r\n', num2str(1/handles.timestep));
        fprintf(fid, '%s\r\n', '%% Total Events Found');
        fprintf(fid, '%s\r\n', get(handles.text_total_number_of_events_found,'String'));
        fprintf(fid, '%s\r\n', '%% END %%');
        fclose(fid); % close file       
    end
    
    fun_make_new_trace(handles, dir_str, 0) % or 0 for normal
    
      
    fun_load_listbox(pwd,handles); % refresh listbox since new directory was created
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                       fun_make_new_trace                          %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a new trace file with all detected events
% Inputs:
%        selectoutput - 0 for nowmal, 1 for OpenNanopore
function fun_make_new_trace(handles, dir_str, selectoutput) % take only the trace for selected events

    indpos = 1; % pointer to current position, initialize
    % if selectoutput == 1
    %     extrapoints = 100;
    % else
        extrapoints = str2num(get(handles.edit_extra_points,'String')); % get number of extra points to each side of the event
    % end
    newtrace = 0;
    newstart = 0;

    num_events_make_new_trace = size(handles.event_info_master,1);

    % pre allocation for speed
    newstart = zeros(num_events_make_new_trace);
    minmax_within_events = zeros(num_events_make_new_trace,10);
    newend = zeros(num_events_make_new_trace);

    for i=1:num_events_make_new_trace % for number of events

        eventpoints = [(handles.event_info_master(i,8)) (handles.event_info_master(i,9))]; % get event start stop

        eventlength = handles.event_info_master(i,9) - handles.event_info_master(i,8) + extrapoints*2 + 1; % calc length in points
        eventtrace = handles.event_trace_master(eventpoints(1)-extrapoints:eventpoints(2)+extrapoints) - handles.event_info_master(i,5); % get the event trace, and remove baseline
        newstart(i) = indpos + extrapoints; % calculate the new start position

        eventtrace_inside = handles.event_trace_master(eventpoints(1):eventpoints(2)) - handles.event_info_master(i,5);

        eventtrace_outside_extra_start = handles.event_trace_master(eventpoints(1)-extrapoints:eventpoints(1)-1) - handles.event_info_master(i,5);
        eventtrace_outside_extra_end = handles.event_trace_master(eventpoints(2)+1:eventpoints(2)+extrapoints) - handles.event_info_master(i,5);

        eventtrace_outside_extra = horzcat(eventtrace_outside_extra_start,eventtrace_outside_extra_end);

        min_in = min(eventtrace_inside);
        max_in = max(eventtrace_inside);
        min_extra = min(eventtrace_outside_extra);
        max_extra = max(eventtrace_outside_extra);
        std_extra = std(eventtrace_outside_extra);
        mean_extra = mean(eventtrace_outside_extra);
        std_start = std(eventtrace_outside_extra_start);
        mean_start = mean(eventtrace_outside_extra_start);
        std_end = std(eventtrace_outside_extra_end);
        mean_end = mean(eventtrace_outside_extra_end);

        minmax_within_events(i,:) = [min_in max_in min_extra max_extra std_extra mean_extra std_start abs(mean_start) std_end abs(mean_end)];

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

        newtrace(indpos:indpos+eventlength-1) = eventtrace; % add this trace to new trace

        newend(i) = newstart(i) + handles.event_info_master(i,9) - handles.event_info_master(i,8); % calculate the new end position
        indpos = indpos + eventlength; % calc new start index value
    end



    if selectoutput == 0
        % write analysis_rawtrace_baseline file
        file_str = strcat(dir_str,'analysis_rawtrace_baseline');
        dlmwrite(file_str, newtrace, 'delimiter', '\n','precision', '%.8f','newline', 'pc');

        % write analysis_events_min_max file
        file_str = strcat(dir_str,'analysis_events_min_max');
        dlmwrite(file_str, minmax_within_events, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');

    end

    if selectoutput == 1
        % write analysis_rawtrace_baseline file
        file_str = strcat(dir_str,'OpenNanopore', filesep, 'events.mat');
        newtrace = transpose(newtrace);
        save(file_str, 'newtrace');

    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_cycle_all_Callback(hObject, eventdata, handles)

    set(handles.button_cycle_all,'BackgroundColor',[1 0 0])

    if get(handles.check_cycle_subfolders,'Value') % cycling through subfolders
        fun_look_for_folders(handles)
    else % analyze just this folder
        fun_cycle_all(handles)
    end
    set(handles.button_cycle_all,'BackgroundColor',[0.925 0.914 0.847])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                            fun_cycle_all                          %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cycles through all files in a folder
function fun_cycle_all(handles) 

    fun_load_listbox(pwd,handles); % load current directory
    listofitems = dir;
    num_files = size(listofitems);
    fun_reset_master_files() % reset everything
    set(handles.text_total_number_of_events_found,'String',num2str(0)) % reset total events counter
    evtotal = 0;
    dirskipped = 0; % # directories skipped
    filecounter = 0; % file counter

    for i= 3:num_files(1,1) % skip first two items, and cycle through rest

        counter_num = i - 2; % same as i but counting from 1
        filezname = char(listofitems(i).name); % get current name

        if listofitems(i).isdir == 0 && filezname(1,1) ~= '~' && ~strcmp(filezname(end-5:end),'_index') % not folder or files starting with ~

            set(handles.listbox_main,'Value',i); % refresh listbox selection to current file

            handles.filename_str = listofitems(i).name;
            filename_str = listofitems(i).name;
            set(handles.edit_status,'String',strcat(handles.filename_str, ' loaded.'))

            fun_load_trace_file(handles.filename_str,handles.text_path,0) % load the file
            handles = guidata(gcbo); % update handles


            if handles.file_load_successful == 1 % file loaded successfully
                type_of_file = get(handles.pop_file_type,'Value');
                if type_of_file == 5 || type_of_file == 6 % if this is an ABF file or TDMS
                    % this is done to be able to cycle over segments

                    filecounter = filecounter + 1; % file is good
                    set(handles.edit_file_number,'String',num2str(filecounter))

                    fun_calculate_resistance()

                    % detect events but dont plot
                    fun_events_detect(0)

                    % update events found counters
                    evfound = str2double(get(handles.text_number_of_events_found,'String'));
                    evtotal = str2double(get(handles.text_total_number_of_events_found,'String'));
                    evtotal = evtotal + evfound;
                    set(handles.text_total_number_of_events_found,'String',num2str(evtotal))

                    drawnow % update GUI

                    % add to master event matrix
                    fun_events_add_to_master_file()

                    % go to second segment
                    handles.currentsegment = 2;

                    % display current segment 
    %                 disp(['Current segment: ' num2str(handles.currentsegment)])


                    while handles.currentsegment < handles.total_number_of_segments + 1 % cycle segments

                        fun_load_trace_file(filename_str,handles.text_path,handles.currentsegment) % load segment

                        handles.currentsegment = get(handles.pop_abf_segment_number,'Value'); % get segment #

                        fun_calculate_resistance()


                        % detect events but dont plot
                        fun_events_detect(0)

                        % update events found counters
                        evfound = str2double(get(handles.text_number_of_events_found,'String'));
                        evtotal = str2double(get(handles.text_total_number_of_events_found,'String'));
                        evtotal = evtotal + evfound;
                        set(handles.text_total_number_of_events_found,'String',num2str(evtotal))

                        drawnow % update GUI

                        % add to master event matrix
                        fun_events_add_to_master_file()

                        % go to next segment
                        handles.currentsegment = handles.currentsegment + 1;

                        % display current segment 
    %                     disp(['Current segment: ' num2str(handles.currentsegment)])

                    end
                    handles.currentsegment = 0; % reset current segment to 0, for next file

                else % every other type of file (Not ABF)

                    filecounter = filecounter + 1; % file is good
                    set(handles.edit_file_number,'String',num2str(filecounter))

                    fun_calculate_resistance()

                    % detect events, do not plot
                    fun_events_detect(0) 
                    drawnow % update GUI

                    % update events found counters
                    evfound = str2double(get(handles.text_number_of_events_found,'String'));
                    evtotal = str2double(get(handles.text_total_number_of_events_found,'String'));
                    evtotal = evtotal + evfound;
                    set(handles.text_total_number_of_events_found,'String',num2str(evtotal))

                    % add to master event matrix
                    fun_events_add_to_master_file()
                end
            end
        else
            dirskipped = dirskipped + 1; % this one got skipped
        end
        fclose('all'); % close stuff
    end


    if get(handles.check_auto_save,'Value') && evtotal > 0 % autosave on and events were found
        handles = guidata(gcbo); % update handles
        fun_save_detection_results(handles) % save it
        fun_reset_master_files() % reset variables
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%             fun_look_for_folders               %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this will build a list of the subfolders in the current path, it will
% then change the path an go inside each subfolder, inside each subfolder
% it will look for a ~detect.par file, if it is there it will load it and
% analyze the trace files within that folder. For each subfolder the
% function calls itself to look for more subfolders, up to an unlimited
% depth.
function fun_look_for_folders(handles)
    folder_counter = 0;
    listofitems = dir;
    num_files = size(listofitems);

    for i= 3:num_files(1,1) % count the number of folders not containing tilda
        foldername = char(listofitems(i).name);
        if listofitems(i).isdir == 1 && foldername(1,1) ~= '~'
            folder_counter = folder_counter +1;
            foldername_array{folder_counter} = listofitems(i).name;
        end
    end
    if folder_counter > 0 % for each folder
        for i= 1:folder_counter
            cd(char(foldername_array(i))) % Jump into folder
            %pwd
            fun_look_for_folders(handles) % recursion
            fun_load_par_file(handles,1) % load parameters, and cycle files

            cd ../ % go back
        end
    end
    fun_load_listbox(pwd,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                fun_load_par_file               %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this loads the ~detect.par file in the local directory and uses it to
% update all of detection parameters
function fun_load_par_file(handles, cycleon)
    filenameparam = strcat(pwd, filesep, '~detect.par');
    if isequal(exist(filenameparam),2) % 
        fid = fopen(filenameparam);
        if fid == -1
            disp(strcat(filenameparam,' could not be loaded!!'))
            %set(handles.button_eventrate,'Enable','off')
        else

            line = fgetl(fid);
            line = fgetl(fid); %% Version #

            line = fgetl(fid);
            set(handles.edit_filter_fc,'String',fgetl(fid)); %% Filter Freq

            line = fgetl(fid);
            set(handles.edit_ab_gain_value,'String',fgetl(fid)); %% Alpha*Beta

            line = fgetl(fid);
            set(handles.edit_voltage,'String',fgetl(fid)); %% Voltage (mV)

            line = fgetl(fid);
            saveunfiltered = str2double(fgetl(fid)); %% Save Unfiltered
            if saveunfiltered
                set(handles.check_save_unfiltered_trace,'Value',1);
            else
                set(handles.check_save_unfiltered_trace,'Value',0);
            end

            line = fgetl(fid);
            saveOpenNanopore = str2double(fgetl(fid)); %% Save OpenNanopore Files
            if saveOpenNanopore
                set(handles.check_make_OpenNanopore_files,'Value',1);
            else
                set(handles.check_make_OpenNanopore_files,'Value',0);
            end

            line = fgetl(fid);
            set(handles.edit_window_size,'String',fgetl(fid)); %% Window Size (points)

            line = fgetl(fid);
            set(handles.edit_peak_detection_factor,'String',fgetl(fid)); %% Peak Det Factor

            line = fgetl(fid);
            set(handles.edit_minimum_dwell_time,'String',fgetl(fid)); %% Min Dwell Time (ms)

            line = fgetl(fid);
            set(handles.edit_maximum_dwell_time,'String',fgetl(fid)); %% Max Dwell Time (ms)

            line = fgetl(fid);
            set(handles.edit_skip_points_between_events,'String',fgetl(fid)); %% Min Events Gap

            line = fgetl(fid);
            set(handles.edit_extra_points,'String',fgetl(fid)); %% Extra Points (points)

            line = fgetl(fid);
            set(handles.edit_number_iterations,'String',fgetl(fid)); %% Number of Iterations

            line = fgetl(fid);
            filttype = str2double(fgetl(fid)); %% Filter Type (1-gauss)
            if filttype ~= 0
                set(handles.pop_filter_type,'Value',filttype);
                set(handles.check_use_filter,'Value',1);
            else
                set(handles.check_use_filter,'Value',0);
            end

            line = fgetl(fid);
            set(handles.edit_sigma_moving_STD_win_size,'String',fgetl(fid)); %% Moving STD Window Size

            line = fgetl(fid);
            filetypevar = str2double(fgetl(fid)); %% File Type
            if filetypevar == 1 || filetypevar == 2 || filetypevar == 3 || filetypevar == 4 || filetypevar == 5 || filetypevar == 6 || filetypevar == 7
                set(handles.pop_file_type,'Value',filetypevar);
            else
                set(handles.pop_file_type,'Value',1);
            end

            line = fgetl(fid);
            baselinemethod = str2double(fgetl(fid)); %% Baseline Method
            if baselinemethod == 1 || baselinemethod == 2 || baselinemethod == 3 || baselinemethod == 4
                set(handles.pop_baseline_method,'Value',baselinemethod);
            else
                set(handles.pop_baseline_method,'Value',1);
            end


            line = fgetl(fid);
            set(handles.edit_manual_baseline_value,'String',fgetl(fid)); %% Manual Baseline Value

            line = fgetl(fid);
            sigmamethod = str2double(fgetl(fid)); %% Sigma Method
            if sigmamethod == 1 || sigmamethod == 2 || sigmamethod == 3 || sigmamethod == 4
                set(handles.pop_sigma_method,'Value',sigmamethod);
            else
                set(handles.pop_sigma_method,'Value',1);
            end

            line = fgetl(fid);
            set(handles.edit_sigma,'String',fgetl(fid)); %% Sigma Value

            line = fgetl(fid);
            set(handles.edit_STD_hist_bins,'String',fgetl(fid)); %% STD Hist Bins

            line = fgetl(fid);
            set(handles.edit_current_hist_bins,'String',fgetl(fid)); %% Current Hist Bins

            line = fgetl(fid);
            event_direction = str2double(fgetl(fid)); %% Event Type (1-down 2-up 3-both)
            if event_direction == 1 || event_direction == 2 || event_direction == 3
                set(handles.pop_event_direction,'Value',event_direction);
            else
                set(handles.pop_event_direction,'Value',1);
            end

        end
        fclose(fid);
        if cycleon == 1
            fun_cycle_all(handles) % detect stuff and save
        end
    else
        set(handles.edit_status,'String','No PAR file found')
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%          fun_events_add_to_master_file         %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add the events detected in the local trace to the master trace
function fun_events_add_to_master_file()
    handles=guidata(gcbo); % get handles

    if isfield(handles,'event_info') % if there are events to add
        listofitems = dir;
        selectnum = get(handles.listbox_main,'Value');
        filename_str = listofitems(selectnum).name; % filename

        if isfield(handles,'event_info_master') % if master file exists already

            mastertracelength = length(handles.event_trace_master);

        else

            mastertracelength = 0;

        end

        event_info_temp = handles.event_info;

        event_info_temp(:,8) = handles.event_start_trace_new + mastertracelength; % update start pointer
        event_info_temp(:,9) = handles.event_stop_trace_new + mastertracelength; % update end pointer

        % TODO, next four lines can be optimized
        event_info_temp(:,10) = handles.event_info(:,8); % pointer of event in original trace
        event_info_temp = sortrows(event_info_temp,10); % sort the events by starting position
        event_info_temp = event_info_temp(:,1:9); % remove column used to sort
        event_info_time = sortrows(handles.event_info,8); % remember how you sorted everything

        handles.event_more_info(:,6) = handles.event_info(:,8); % add the starting points to more_info matrix so it can be sorted 
        handles.event_more_info = sortrows(handles.event_more_info,6); % sort it the same way you sorted the last one
        handles.event_more_info = handles.event_more_info(:,1:5); % remove the starting points column

        extra_data_fields = 4; % add four more columns
        temp_storage = cell(length(handles.event_start_trace_new),extra_data_fields);

        if get(handles.pop_file_type,'Value') == 5 % Axon binaries
            temp_storage(:,1) = {strcat('"',filename_str,' segment ',get(handles.edit_segment_number,'String'),'"')}; % name of the file, append segment #
            % file number*1000000 + segment number
            fileidnum = 1000000*str2double(get(handles.edit_file_number,'String')) + str2double(get(handles.edit_segment_number,'String'));
            temp_storage(:,2) = cellstr(num2str(fileidnum));
        else % everything else
            temp_storage(:,1) = {strcat('"',filename_str,'"')}; % name of the file
            temp_storage(:,2) = cellstr(num2str(get(handles.edit_file_number,'String'))); % file number
        end

        temp_storage(:,3) = {get(handles.edit_resistance,'String')}; % resistance
        temp_storage(:,4) = cellstr(num2str(handles.timestep*event_info_time(:,8))); % time local
        totalfiletime = length(handles.filtered_trace)*handles.timestep;

        if ~isfield(handles,'total_time_so_far')
            handles.total_time_so_far = 0;
        end

        temp_storage(:,5) = cellstr(num2str(handles.timestep*event_info_time(:,8)+handles.total_time_so_far)); % time total

        % inside temp storage:
        % 1 name of file
        % 2 file number
        % 3 resistance
        % 4 time - local
        % 5 time - global

        handles.total_time_so_far = handles.total_time_so_far + totalfiletime; % add length of this file, for next time

        if isfield(handles,'event_info_master') % if master file exists already

            handles.event_info_master = vertcat(handles.event_info_master,event_info_temp); % add the new event info to the previous ones in master trace
            handles.event_trace_master = [handles.event_trace_master handles.events_local_trace]; % add local trace to master
            if get(handles.check_save_unfiltered_trace,'Value') % save unfiltered trace
                handles.event_unfiltered_trace_master = [handles.event_unfiltered_trace_master handles.event_unfiltered_trace];
            end
            handles.event_more_info_master = vertcat(handles.event_more_info_master,handles.event_more_info); % add the new more_info
            handles.event_master_trace_info = vertcat(handles.event_master_trace_info,temp_storage); % add the trace info
        else

            handles.event_info_master = event_info_temp; % create it
            handles.event_trace_master = handles.events_local_trace; % create it
            if get(handles.check_save_unfiltered_trace,'Value') % save unfiltered trace
                handles.event_unfiltered_trace_master = handles.event_unfiltered_trace;
            end
            handles.event_more_info_master = handles.event_more_info; % create it
            handles.event_master_trace_info = temp_storage; % create it
        end

        set(handles.edit_status,'String',strcat(num2str(length(handles.event_info_master)), ' events in the master file.'))
    end
    guidata(gcbo, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_find_events_Callback(hObject, eventdata, handles)
    % find events
    fun_events_detect(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   button_add_local_events_to_master_trace_Callback    %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manually add the detected local trace events to the master trace
function button_add_local_events_to_master_trace_Callback(hObject, eventdata, handles)

    evfound = str2double(get(handles.text_number_of_events_found,'String'));
    evtotal = str2double(get(handles.text_total_number_of_events_found,'String'));
    evtotal = evtotal + evfound;
    set(handles.text_total_number_of_events_found,'String',num2str(evtotal)) % update the master events total
    fun_events_add_to_master_file()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_reset_master_event_trace_Callback(hObject, eventdata, handles)
    % erase all of the master variables
    fun_reset_master_files()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%               fun_reset_master_files           %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reset all of the master file variables
function fun_reset_master_files
    handles=guidata(gcbo);
    if isfield(handles,'event_info_master')
        handles = rmfield (handles,'event_info_master');
        handles = rmfield (handles,'event_master_trace_info');
    end
    if isfield(handles,'event_trace_master')
        handles = rmfield (handles,'event_trace_master');
    end
    if isfield(handles,'event_unfiltered_trace_master')
        handles = rmfield (handles,'event_unfiltered_trace_master');
    end
    if isfield(handles,'event_more_info_master')
        handles = rmfield (handles,'event_more_info_master');
    end
    if isfield(handles,'total_time_so_far')
        handles = rmfield (handles,'total_time_so_far');
    end
    set(handles.text_total_number_of_events_found,'String',num2str(0))
    set(handles.edit_segment_number,'String',num2str(0))
    guidata(gcbo, handles); % update handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text_number_of_events_found_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text_number_of_events_found_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text_total_number_of_events_found_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function text_total_number_of_events_found_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider_events_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider_events_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%       edit_event_number_event_info_Callback    %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot some event if its number is entered into the text box
function edit_event_number_event_info_Callback(hObject, eventdata, handles)
    % used to view different events, in the local browser
    numevent = str2num(get(handles.edit_event_number_event_info,'String'));
    % if number makes sense
    if get(handles.radio_browse_local_events,'Value') == 1 % this is for browsing local events
        totalevents = str2num(get(handles.text_number_of_events_found,'String'));
        if (numevent > totalevents) || (numevent < 1) % sanity check
            numevent = 1; 
        end
    else % this is for browsing master events
        totalevents = str2num(get(handles.text_total_number_of_events_found,'String'));
        if (numevent > totalevents) || (numevent < 1) % sanity check
            numevent = 1;
        end
    end
    set(handles.edit_event_number_event_info,'String',num2str(numevent));

    if totalevents ~= 0
        fun_event_plot('here') % plot it
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_event_number_event_info_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_amplitude_event_info_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_amplitude_event_info_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_integral_event_info_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_integral_event_info_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_max_amplitude_event_info_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_max_amplitude_event_info_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_dwell_event_info_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_dwell_event_info_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_baseline_event_info_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_baseline_event_info_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_detection_level_event_info_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_detection_level_event_info_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_start_points_event_info_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_start_points_event_info_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_stop_points_event_info_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_stop_points_event_info_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_open_Callback(hObject, eventdata, handles)
    % open new path
    pathvar = uigetdir;
    cd (pathvar);
    fun_load_listbox(pathvar,handles)
    guidata(hObject, handles); % update handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_plot_trace_in_win_Callback(hObject, eventdata, handles)
    % open trace plot in new window
    fun_trace_plot('new',0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_plot_event_Callback(hObject, eventdata, handles)
    % open event plot in new window
    fun_event_plot('new')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_file_type_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_file_type_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_status_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_status_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_diagnostic_plot_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_plotevents_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_zoom_Callback(hObject, eventdata, handles)
    % start  case zoom 1 pan 0
    % switch case zoom 0 pan 1
    handles.zoomstring = get(handles.button_zoom,'String');
    handles.panstring = get(handles.button_pan,'String');

    if strcmp(handles.zoomstring, 'Zoom On')
        % we are in Zoom On mode and want to switch
        zoom off
        set(handles.button_zoom,'String','Zoom Off');
    else
        % we are in Zoom Off mode and want to switch
        pan off
        set(handles.button_pan,'String','Pan Off');
        zoom on
        set(handles.button_zoom,'String','Zoom On');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pan_Callback(hObject, eventdata, handles)

    handles.zoomstring = get(handles.button_zoom,'String');
    handles.panstring = get(handles.button_pan,'String');

    if strcmp(handles.panstring, 'Pan On')
        % we are in Pan On mode and want to switch
        pan off
        set(handles.button_pan,'String','Pan Off');
    else
        % we are in Pan Off mode and want to switch
        zoom off
        set(handles.button_zoom,'String','Zoom Off');
        pan on
        set(handles.button_pan,'String','Pan On');
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_event_browse_up_event_info_Callback(hObject, eventdata, handles)
    % browse through the different events
    numevent = str2num(get(handles.edit_event_number_event_info,'String'));
    if get(handles.radio_browse_local_events,'Value') == 1 % this is for browsing local events
        totalevents = str2num(get(handles.text_number_of_events_found,'String'));
        if numevent < totalevents  % sanity check
            numevent = numevent + 1;
        end
    else  % this is for browsing master events
        totalevents = str2num(get(handles.text_total_number_of_events_found,'String'));
        if numevent < totalevents  % sanity check
            numevent = numevent + 1;
        end
    end
    set(handles.edit_event_number_event_info,'String',num2str(numevent));
    fun_event_plot('here') % plot it

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_event_browse_down_event_info_Callback(hObject, eventdata, handles)
    % browse through the different events
    numevent = str2num(get(handles.edit_event_number_event_info,'String'));
    if numevent > 1 % sanity check
        numevent = numevent - 1;
    end
    if get(handles.radio_browse_local_events,'Value') == 1
        totalevents = str2num(get(handles.text_number_of_events_found,'String'));
        if numevent > totalevents % sanity check
            numevent = totalevents;
        end
    else
        totalevents = str2num(get(handles.text_total_number_of_events_found,'String'));
        if numevent > totalevents % sanity check
            numevent = totalevents;
        end
    end
    set(handles.edit_event_number_event_info,'String',num2str(numevent));
    fun_event_plot('here') % plot it

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                   fun_event_plot               %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this will plot a single detected event in the axis in the top right
% Inputs:
% location - 'here' to plot inside GUI, something else to plot in a new
%                   figure
function fun_event_plot(location)
    handles=guidata(gcbo);
    numevent = str2num(get(handles.edit_event_number_event_info,'String')); % what event are we looking at
    extradisplay = str2num(get(handles.edit_display_extra_points_event_info,'String')); % how many extra points to show
    
    
    if get(handles.radio_browse_local_events,'Value') == 1 % local events
        
        startpoint_var = handles.event_info(numevent,8);
        stoppoint_var = handles.event_info(numevent,9);
        
        amplitude_var = handles.event_info(numevent,1);
        
        q2_var = handles.event_more_info(numevent,2);
        q3_var = handles.event_more_info(numevent,3);
    
        eventpoints = [(startpoint_var) (stoppoint_var)]; % the start and stop points
        
        %create a time vector
        timetemp = transpose(0:(handles.timestep*1000):handles.timestep*(eventpoints(2)-eventpoints(1)+2*extradisplay)*1000);
        
        startstoptemp_time = [(handles.timestep*1000*(extradisplay)) handles.timestep*1000*((eventpoints(2)-eventpoints(1))+extradisplay)];
        
        eventtrace = handles.filtered_trace(eventpoints(1)-extradisplay:eventpoints(2)+extradisplay);
        startstoptemp_value = [handles.filtered_trace(eventpoints(1)) handles.filtered_trace(eventpoints(2))];
        temprunavg = handles.trace_runavg_fwd(eventpoints(1)).*ones(size(timetemp));
        
        if handles.event_info(numevent,7) == 1 % up events
            tempamplitude = (handles.trace_runavg_fwd(eventpoints(1)) + amplitude_var).*ones(size(timetemp));
        else % down events
            tempamplitude = (handles.trace_runavg_fwd(eventpoints(1)) - amplitude_var).*ones(size(timetemp));
        end
        if (get(handles.check_remove_baseline,'Value') == 1) ||  (get(handles.check_conductance,'Value') == 1)% baseline is shifted to 0 nA
            eventtrace = eventtrace - handles.trace_runavg_fwd(eventpoints(1));
            startstoptemp_value = startstoptemp_value - handles.trace_runavg_fwd(eventpoints(1)); 
            temprunavg = temprunavg - handles.trace_runavg_fwd(eventpoints(1));
            tempamplitude = tempamplitude - handles.trace_runavg_fwd(eventpoints(1));
        end

        

        set(handles.edit_max_amplitude_event_info,'String',num2str(handles.event_info(numevent,2)));
        set(handles.edit_integral_event_info,'String',num2str(handles.event_info(numevent,3)));
        set(handles.edit_dwell_event_info,'String',num2str(handles.event_info(numevent,4)));
        set(handles.edit_baseline_event_info,'String',num2str(handles.event_info(numevent,5)));
        set(handles.edit_detection_level_event_info,'String',num2str(handles.event_info(numevent,6)));

    % 1average	2maximum	3integral	4dwelltime	5baseline	6detectionlevel
                        % 7eventtype	8start point	9end point
    else %%%%%%%%%%%% Event in Master Trace
        
        startpoint_var = handles.event_info_master(numevent,8);
        stoppoint_var = handles.event_info_master(numevent,9);
        amplitude_var = handles.event_info_master(numevent,1);
        q2_var = handles.event_more_info_master(numevent,2);
        q3_var = handles.event_more_info_master(numevent,3);
        baseline_var = handles.event_info_master(numevent,5);
        
        
        extrapointsintrace = str2num(get(handles.edit_extra_points,'String'));
        if extradisplay > extrapointsintrace % cant show more points than we have in the trace
            extradisplay = extrapointsintrace;
        end
        % the start an stop points
        eventpoints = [(startpoint_var+1) (stoppoint_var+1)]; % TODO why is this plus one?
        
        %create a time vector
        timetemp = transpose(0:(handles.timestep*1000):handles.timestep*(eventpoints(2)-eventpoints(1)+2*extradisplay)*1000);
        
        % starting and stopping point times
        startstoptemp_time = [(handles.timestep*1000*(extradisplay)) handles.timestep*1000*((eventpoints(2)-eventpoints(1))+extradisplay)];
        
        % get the event points
        eventtrace = handles.event_trace_master(eventpoints(1)-extradisplay:eventpoints(2)+extradisplay);
        
        % values of the starting and stoping points
        startstoptemp_value = [handles.event_trace_master(eventpoints(1)) handles.event_trace_master(eventpoints(2))];
        
        % the baseline
        temprunavg = baseline_var.*ones(size(eventtrace));
        
        % amplitude (positive for down, neg for up)
        tempamplitude = -1.*((baseline_var - amplitude_var).*ones(size(eventtrace)));
        
        if (get(handles.check_remove_baseline,'Value') == 1) ||  (get(handles.check_conductance,'Value') == 1)% baseline is shifted to 0 nA
            eventtrace = eventtrace - temprunavg(1);
            startstoptemp_value = startstoptemp_value - temprunavg(1);
            temprunavg = temprunavg - temprunavg(1);
            
            % this is the line showing amplitude of the event
            if handles.event_info_master(numevent,7) == 1 % up events
                tempamplitude = (0 + amplitude_var)*ones(size(eventtrace));
            else % down events
                tempamplitude = (0 - amplitude_var)*ones(size(eventtrace));
            end
        end
        
        set(handles.edit_max_amplitude_event_info,'String',num2str(handles.event_info_master(numevent,2)));
        set(handles.edit_integral_event_info,'String',num2str(handles.event_info_master(numevent,3)));
        set(handles.edit_dwell_event_info,'String',num2str(handles.event_info_master(numevent,4)));
        set(handles.edit_detection_level_event_info,'String',num2str(handles.event_info_master(numevent,6)));
        set(handles.edit_baseline_event_info,'String',num2str(baseline_var));

    % 1average	2maximum	3integral	4dwelltime	5baseline	6detectionlevel
                        % 7eventtype	8start point	9end point
    end
    
    if get(handles.check_conductance,'Value') == 1 % trace converted to nS
        voltage = str2double(get(handles.edit_voltage,'String'))/1000; % in Volts
        eventtrace = eventtrace./voltage;
        startstoptemp_value = startstoptemp_value./voltage;  
        temprunavg = temprunavg./voltage;
        tempamplitude = tempamplitude./voltage;
    end
    
    if strcmp(location, 'here')
        cla(handles.axes_event,'reset')
        axes(handles.axes_event)
    else
        figure()
    end
    
    if get(handles.radio_browse_local_events,'Value') == 1 % local events
        
        if get(handles.check_diagnostic_plot,'Value') == 1
            plot(timetemp, temprunavg, '-m'); % running average line
            hold on;
            plot(timetemp, tempamplitude, '-r'); % amplitude line      
            plot(timetemp(q2_var+extradisplay), eventtrace(q2_var+extradisplay), '.b','MarkerSize',12);
            plot(timetemp(q3_var+extradisplay), eventtrace(q3_var+extradisplay), '.r','MarkerSize',12);
            hold on;
        end
    else
        
        if get(handles.check_diagnostic_plot,'Value') == 1
            plot(timetemp, temprunavg, '-m'); % running average line
            hold on;
            plot(timetemp, tempamplitude, '-r'); % amplitude line
            plot(timetemp(q2_var+extradisplay-1), eventtrace(q2_var+extradisplay-1), '.b','MarkerSize',12); % TODO why -1
            plot(timetemp(q3_var+extradisplay-1), eventtrace(q3_var+extradisplay-1), '.r','MarkerSize',12); % TODO why -1
            hold on;
        end
        
    end
    
    set(handles.edit_amplitude_event_info,'String',num2str(amplitude_var));
    set(handles.edit_start_points_event_info,'String',num2str(startpoint_var));
    set(handles.edit_stop_points_event_info,'String',num2str(stoppoint_var));
        
    set(handles.edit_start_time_event_info,'String',num2str(handles.timestep*startpoint_var));
    set(handles.edit_stop_time_event_info,'String',num2str(handles.timestep*stoppoint_var));
    
    plot(startstoptemp_time, startstoptemp_value, '.k','MarkerSize',20);
    hold on;
    plot(timetemp,eventtrace,'-b')
    
    ylimits = get(handles.axes_event,'YLim');
    xlimits = get(handles.axes_event,'XLim');
    axis([0 max(timetemp) ylimits(1) ylimits(2)])

    hXLabel = xlabel('Time (ms)');
    if get(handles.check_conductance,'Value') == 0
        hYLabel = ylabel('Current (nA)');
    else
        hYLabel = ylabel('Conductance (nS)');
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
      'XColor'      , [0 0 0], ...
      'YColor'      , [0 0 0], ...
      'LineWidth'   , 1         );
    
    guidata(gcbo, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%            fun_calculate_resistance            %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the resistance of the currently loaded trace in Mohms
function fun_calculate_resistance()
    handles=guidata(gcbo);
    voltage = str2num(get(handles.edit_voltage,'String'))/1000; % in Volts
    handles.resistance = (voltage/(handles.trace_mean*1E-9))/1E6; % resistance in Mohms
    set(handles.edit_resistance,'String',num2str(handles.resistance)); % update textbox
    guidata(gcbo, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_display_extra_points_event_info_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_display_extra_points_event_info_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_start_time_event_info_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_start_time_event_info_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_stop_time_event_info_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_stop_time_event_info_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_highlight_events_Callback(hObject, eventdata, handles)
    % empty



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_remove_baseline_Callback(hObject, eventdata, handles)
    % empty



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_conductance_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_voltage_Callback(hObject, eventdata, handles)
    % update the resistance value
    fun_calculate_resistance()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_voltage_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_resistance_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_resistance_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_plot_trace_in_gui_Callback(hObject, eventdata, handles)
    % update trace
    fun_trace_plot('here',0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%    button_dump_variables_to_matlab_Callback    %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this will dump variables to Matlab base workspace
function button_dump_variables_to_matlab_Callback(hObject, eventdata, handles)

    if isfield(handles,'event_info')
        assignin('base', 'event_info', handles.event_info)
    end
    if isfield(handles,'sigma_window')
        assignin('base', 'sigma_window', handles.sigma_window)
    end
    if isfield(handles,'event_master_trace_info')
        assignin('base', 'event_master_trace_info', handles.event_master_trace_info)
    end
    if isfield(handles,'trace')
        assignin('base', 'trace', handles.trace)
    end
    if isfield(handles,'time_vector')
        assignin('base', 'time_vector', handles.time_vector)
    end
    if isfield(handles,'cutoff_frequency')
        assignin('base', 'cutoff_frequency', handles.cutoff_frequency)
    end
    if isfield(handles,'filtered_trace')
        assignin('base', 'filtered_trace', handles.filtered_trace)
    end
    if isfield(handles,'trace_runavg_fwd')
        assignin('base', 'trace_runavg_fwd', handles.trace_runavg_fwd)
    end
    if isfield(handles,'resistance')
        assignin('base', 'resistance', handles.resistance)
    end
    if isfield(handles,'event_trace_master')
        assignin('base', 'event_trace_master', handles.event_trace_master)
    end
    if isfield(handles,'event_info_master')
        assignin('base', 'event_info_master', handles.event_info_master)
    end
    if isfield(handles,'event_stop_trace_new')
        assignin('base', 'event_stop_trace_new', handles.event_stop_trace_new)
    end
    if isfield(handles,'event_start_trace_new')
        assignin('base', 'event_start_trace_new', handles.event_start_trace_new)
    end
    if isfield(handles,'event_time_withoutextra')
        assignin('base', 'event_time_withoutextra', handles.event_time_withoutextra)
    end
    if isfield(handles,'events_local_trace_withoutextra')
        assignin('base', 'events_local_trace_withoutextra', handles.events_local_trace_withoutextra)
    end
    if isfield(handles,'event_startstop_time')
        assignin('base', 'event_startstop_time', handles.event_startstop_time)
    end
    if isfield(handles,'event_startstop_trace')
        assignin('base', 'event_startstop_trace', handles.event_startstop_trace)
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_file_number_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_file_number_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_make_OpenNanopore_files_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_master_in_win_Callback(hObject, eventdata, handles)
    % plot master trace in a new window
    fun_master_file_plot('new')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_master_in_gui_Callback(hObject, eventdata, handles)
    % plot master trace here
    fun_master_file_plot('here')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%               fun_master_file_plot             %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this will plot all concatenated events in the master file
function fun_master_file_plot(location) % plotting the master trace
    handles=guidata(gcbo); % get handles
    
    if strcmp(location, 'here') % plot it here in the GUI
        cla(handles.axes_trace,'reset')
        axes(handles.axes_trace)
    end
    
    if strcmp(location, 'new') % plot it in a new figure
        figure()
    end
    
    if isfield(handles,'event_trace_master') % the master trace exists
        
        lengthtrace = length(handles.event_trace_master);
        temptime = 0:handles.timestep:(lengthtrace-1)*handles.timestep; % create the time vector
        
        startvalue = handles.event_trace_master(handles.event_info_master(:,8)); % all event start values
        starttime = temptime(handles.event_info_master(:,8)); % all event start times
        
        endvalue = handles.event_trace_master(handles.event_info_master(:,9)); % all event end values
        endtime = temptime(handles.event_info_master(:,9)); % all event end times
        
        plot(temptime,handles.event_trace_master,'b-') % show the master trace
        hold on
        
        plot(starttime, startvalue, '.k','MarkerSize',5); % plot start points
        plot(endtime, endvalue, '.r','MarkerSize',5); % plot end points
        
        axis([0 max(temptime) min(handles.event_trace_master) max(handles.event_trace_master)]) % set limits
    end
    guidata(gcbo, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_browse_local_events_Callback(hObject, eventdata, handles)

    % Hint: get(hObject,'Value') returns toggle state of radio_browse_local_events
    if get(handles.radio_browse_local_events,'Value') % local on
        set(handles.radio_browse_master_events,'Value',0) % master off
    end
    if str2double(get(handles.text_total_number_of_events_found,'String')) == 0
        set(handles.radio_browse_master_events,'Value',0) % master off
        set(handles.radio_browse_local_events,'Value',1) % local on
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_browse_master_events_Callback(hObject, eventdata, handles)

    % Hint: get(hObject,'Value') returns toggle state of radio_browse_master_events
    if get(handles.radio_browse_master_events,'Value') % master on
        set(handles.radio_browse_local_events,'Value',0) % local off
    end
    if str2double(get(handles.text_total_number_of_events_found,'String')) == 0
        set(handles.radio_browse_master_events,'Value',0) % master off
        set(handles.radio_browse_local_events,'Value',1) % local on
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_timestep_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_timestep_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_save_unfiltered_trace_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%              listbox_main_Callback             %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this will move a selected file into a subdirectory named ~unused_traces
function button_move_file_to_unused_Callback(hObject, eventdata, handles)
        
        index_selected = get(handles.listbox_main,'Value');
        file_list = get(handles.listbox_main,'String');
        % Item selected in list box
        filename = file_list{index_selected};
        
        if  isdir(filename) % If folder
            % do nothing
        else
            [pathva,name,ext] = fileparts(filename);
            dir_str = strcat(pwd, filesep, '~unused_traces');
            if(~exist(dir_str, 'dir')) % create ~unused_traces if it does not exist
                mkdir(dir_str);
            end
            fclose('all');
            newname = strcat(dir_str, filesep, filename);
            movefile(filename,newname); % move the file
            
        end
        fun_load_listbox(pwd,handles); % Reload list box
        guidata(gcbo, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_auto_save_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_cycle_subfolders_Callback(hObject, eventdata, handles)

    if get(handles.check_cycle_subfolders,'Value') % subfolder on
        set(handles.check_auto_save,'Value',1) % save on
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_parfile_create_Callback(hObject, eventdata, handles)

    dir_str = strcat(pwd, filesep);
    fun_save_par_file(handles, dir_str)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                 fun_save_par_file              %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this saves a file named ~detect.par in the current path which contains
% all of the detection parameters as they appear in the GUI
function fun_save_par_file(handles, dir_str)
    file_str = strcat(dir_str,'~detect.par');
    fid = fopen(file_str, 'w');
    fprintf(fid, '%s\r\n', '%% Detect Par File');
    fprintf(fid, '%s\r\n', strcat('%% Version ', get(handles.text_version,'String')));
    fprintf(fid, '%s\r\n', '%% Filter Freq');
    fprintf(fid, '%s\r\n', get(handles.edit_filter_fc,'String'));
    fprintf(fid, '%s\r\n', '%% Alpha*Beta');
    fprintf(fid, '%s\r\n', get(handles.edit_ab_gain_value,'String'));
    fprintf(fid, '%s\r\n', '%% Voltage (mV)');
    fprintf(fid, '%s\r\n', get(handles.edit_voltage,'String'));
    fprintf(fid, '%s\r\n', '%% Save Unfiltered');
    fprintf(fid, '%s\r\n', num2str(get(handles.check_save_unfiltered_trace,'Value')));
    fprintf(fid, '%s\r\n', '%% Save OpenNanopore Files');
    fprintf(fid, '%s\r\n', num2str(get(handles.check_make_OpenNanopore_files,'Value')));
    fprintf(fid, '%s\r\n', '%% Window Size (points)');
    fprintf(fid, '%s\r\n', get(handles.edit_window_size,'String'));
    fprintf(fid, '%s\r\n', '%% Peak Det Factor');
    fprintf(fid, '%s\r\n', get(handles.edit_peak_detection_factor,'String'));
    fprintf(fid, '%s\r\n', '%% Min Dwell Time (ms)');
    fprintf(fid, '%s\r\n', get(handles.edit_minimum_dwell_time,'String'));
    fprintf(fid, '%s\r\n', '%% Max Dwell Time (ms)');
    fprintf(fid, '%s\r\n', get(handles.edit_maximum_dwell_time,'String'));
    fprintf(fid, '%s\r\n', '%% Min Events Gap');
    fprintf(fid, '%s\r\n', get(handles.edit_skip_points_between_events,'String'));
    fprintf(fid, '%s\r\n', '%% Extra Points (points)');
    fprintf(fid, '%s\r\n', get(handles.edit_extra_points,'String'));
    fprintf(fid, '%s\r\n', '%% Number of Iterations');
    fprintf(fid, '%s\r\n', get(handles.edit_number_iterations,'String'));
    fprintf(fid, '%s\r\n', '%% Filter Type (0-none 1-gauss 2-butter2 3-butter4)');
    if get(handles.check_use_filter,'Value') == 1
        fprintf(fid, '%s\r\n', num2str(get(handles.pop_filter_type,'Value')));
    else
        fprintf(fid, '%s\r\n', '0');
    end
    fprintf(fid, '%s\r\n', '%% Moving STD Win Size');
    fprintf(fid, '%s\r\n', get(handles.edit_sigma_moving_STD_win_size,'String'));
    fprintf(fid, '%s\r\n', '%% File Type');
    fprintf(fid, '%s\r\n', num2str(get(handles.pop_file_type,'Value')));
    fprintf(fid, '%s\r\n', '%% Baseline Method');
    fprintf(fid, '%s\r\n', num2str(get(handles.pop_baseline_method,'Value')))
    fprintf(fid, '%s\r\n', '%% Baseline Value');
    fprintf(fid, '%s\r\n', get(handles.edit_manual_baseline_value,'String'));
    fprintf(fid, '%s\r\n', '%% Sigma Method');
    fprintf(fid, '%s\r\n', num2str(get(handles.pop_sigma_method,'Value')));
    fprintf(fid, '%s\r\n', '%% Sigma Value');
    fprintf(fid, '%s\r\n', get(handles.edit_sigma,'String'));
    fprintf(fid, '%s\r\n', '%% STD Hist Bins');
    fprintf(fid, '%s\r\n', get(handles.edit_STD_hist_bins,'String'));
    fprintf(fid, '%s\r\n', '%% Current Hist Bins');
    fprintf(fid, '%s\r\n', get(handles.edit_current_hist_bins,'String'));
    fprintf(fid, '%s\r\n', '%% Event Type (1-down 2-up 3-both)');
    fprintf(fid, '%s\r\n', num2str(get(handles.pop_event_direction,'Value')));
    fprintf(fid, '%s\r\n', '%% END PAR FILE %%');
    fclose(fid);
    fun_load_listbox(pwd, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_parfile_load_Callback(hObject, eventdata, handles)

    fun_load_par_file(handles, 0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in pop_analysis_type.
function pop_analysis_type_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_analysis_type_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in pop_filter_type.
function pop_filter_type_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_filter_type_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_stdwin_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_stdwin_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_std_win_size_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_std_win_size_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_manual_baseline_value_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_manual_baseline_value_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%        pop_abf_segment_number_Callback         %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is called when the segment number is changed in the drop down menu.
% The selected current segment is opened.
function pop_abf_segment_number_Callback(hObject, eventdata, handles)

    index_selected = get(handles.listbox_main,'Value');
    file_list = get(handles.listbox_main,'String');

    % Item selected in list box
    filename = file_list{index_selected};
    currentsegment_old = handles.currentsegment;
    handles.currentsegment = get(handles.pop_abf_segment_number,'Value'); % what segment are we on

    if handles.currentsegment ~= currentsegment_old % a new segment is selected
        fun_load_trace_file(filename, handles.text_path,handles.currentsegment) % load the segment
        handles=guidata(gcbo);
        if handles.file_load_successful == 1
            fun_calculate_resistance() % update the resistance
            set(handles.check_highlight_events,'Value',0);
            fun_trace_plot('here',0) % plot it
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_abf_segment_number_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%         button_noise_analysis_Callback         %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate power spectrum for current trace. Uses unfiltered data. See GUI_noise
function button_noise_analysis_Callback(~, ~, handles)
    % plot FFT of current trace

    FsG = 1/(handles.timestep);
    [fGn PxxGn] = calcPSD(FsG, handles.trace);

    noise_fig_handles = figure();
    loglog(fGn,PxxGn,'-r')
    xlabel('Frequency (Hz)')
    ylabel('PSD (nA^2/Hz)') 
    
%   axis([1 3000000 1E-18 1E-4])
    
%    xlimits = get(noise_fig_handles,'XLim');
%    ylimits = get(noise_fig_handles,'YLim');
% 
%    axis([0 5E6 ylimits(1) ylimits(2)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_segment_number_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_segment_number_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GUI_detect_CloseRequestFcn(hObject, eventdata, handles)

    % clean up shared data
    HandleMainGUI=getappdata(0,'HandleMainGUI');

    % Remember that your GUIs might try to getappdata that doesn't exist, you should first test if it does exist
    if (isappdata(0,'HandleMainGUI') && isappdata(HandleMainGUI,'MySharedData'))
        rmappdata(HandleMainGUI,'MySharedData') % do rmappdata for all data shared 
    else
    %do something else, maybe loading default values into those variables
    end

    % Hint: delete(hObject) closes the figure
    delete(hObject);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%            button_clampexiv_Callback           %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generates an IV plot and data for a Clampex IV file.
%
function button_clampexiv_Callback(hObject, eventdata, handles)

        index_selected = get(handles.listbox_main,'Value');
        file_list = get(handles.listbox_main,'String');
        % Item selected in list box
        filename = file_list{index_selected};
        
        if  isdir(filename) % If folder
            % do nothing
        else % assume it is a calmpex IV file
            [pathva,name,ext] = fileparts(filename);
            
            % make new subdir to store data
            dir_str = strcat(pwd, filesep,'~IV');
            if(~exist(dir_str, 'dir'))
                mkdir(dir_str);
            end
            fclose('all');
            
            % call IV script
            [current_data, voltage_data, numsweep, numdp]= plot_IV(filename,0);
            
            % IV plot
            newFig = figure;
            plot(voltage_data, current_data./1000,'.-r')
            xlabel('Voltage (mV)')
            ylabel('Current (nA)')
            saveas(newFig, strcat(pwd, filesep, '~IV', filesep, name, '_plot_IV.fig'), 'fig'); % save as .fig
            saveas(newFig, strcat(pwd, filesep, '~IV', filesep, name, '_plot_IV_large.png'), 'png');
            set(newFig, 'PaperPositionMode', 'auto');
            saveas(newFig, strcat(pwd, filesep, '~IV', filesep, name, '_plot_IV_small.png'), 'png');
            close(newFig)
            
            % save IV data
            file_str = strcat(dir_str, filesep, name, '_IV.dat');
            dlmwrite(file_str, [voltage_data current_data], 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
            
        end
        fun_load_listbox(pwd,handles); % Load list box, since new dir was made
        guidata(gcbo, handles);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_sigma_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_sigma_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_segment_size_in_sec_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_segment_size_in_sec_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%        button_current_hist_make_Callback       %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this generates a histogram of all current points in the trace segment
% currently loaded. This is best for stable baselines. Unstable baselines
% will lead to broad peaks.
function button_current_hist_make_Callback(hObject, eventdata, handles)

    bins_in_current_hist = str2double(get(handles.edit_current_hist_bins,'String'));
    figure()
    
    [n,xout] = hist(handles.filtered_trace(handles.windowSize:end),bins_in_current_hist); % make the histogram
    hist(handles.filtered_trace(handles.windowSize:end),bins_in_current_hist)
        
    % next four lines are needed to make the histogram show properly under
    % log scale
    ph = get(gca,'children');
    vn = get(ph,'Vertices');
    vn(:,2) = vn(:,2) + 1;
    set(ph,'Vertices',vn);
    
    set(gca,'YScale','log') % use log scale on Y
    hXLabel = xlabel('Current (nA)');
    
    hYLabel = ylabel('Counts'); % Y label
    
    % style
    set( gca                       , ...
        'FontName'   , 'Helvetica' );
    set([hXLabel, hYLabel], ...
        'FontName'   , 'AvantGarde');
    set([hXLabel, hYLabel]  , ...
        'FontSize'   , 10          );

    set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'on'      , ...
      'YMinorTick'  , 'on'      , ...
      'YGrid'       , 'off'      , ...
      'XColor'      , [0 0 0], ...
      'YColor'      , [0 0 0], ...
      'LineWidth'   , 1         );




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_number_iterations_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_number_iterations_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_baseline_method_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_baseline_method_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_sigma_method_Callback(hObject, eventdata, handles)
    % empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_sigma_method_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_current_hist_bins_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_current_hist_bins_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_sigma_moving_STD_win_size_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_sigma_moving_STD_win_size_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%          button_STD_histogram_Callback         %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this button generates a figure with a histogram of STD values calculated
% for a moving window across the trace
%
function button_STD_histogram_Callback(hObject, eventdata, handles)


    bins_in_STD_hist = str2double(get(handles.edit_STD_hist_bins,'String'));
    figure()
    
    moving_STD_window_size = 1000; % keep it small
    sigma_moving = movingstd(handles.filtered_trace(400:end),moving_STD_window_size,'central');
    
    [n,xout] = hist(sigma_moving,bins_in_STD_hist); % make the histogram
    hist(sigma_moving,bins_in_STD_hist)
        
    % next four lines are needed to make the histogram show properly under
    % log scale
    ph = get(gca,'children');
    vn = get(ph,'Vertices');
    vn(:,2) = vn(:,2) + 1;
    set(ph,'Vertices',vn);
    
    set(gca,'YScale','log') % use log scale on Y
    hXLabel = xlabel('Sigma (nA)');
    
    hYLabel = ylabel('Counts'); % Y label
    
    % style
    set( gca                       , ...
        'FontName'   , 'Helvetica' );
    set([hXLabel, hYLabel], ...
        'FontName'   , 'AvantGarde');
    set([hXLabel, hYLabel]  , ...
        'FontSize'   , 10          );

    set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'on'      , ...
      'YMinorTick'  , 'on'      , ...
      'YGrid'       , 'off'      , ...
      'XColor'      , [0 0 0], ...
      'YColor'      , [0 0 0], ...
      'LineWidth'   , 1         );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_STD_hist_bins_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_STD_hist_bins_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_event_direction_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pop_event_direction_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_downsample_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_downsample_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_downsample_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_target_res_nA_Callback(hObject, eventdata, handles)
    % empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_target_res_nA_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
