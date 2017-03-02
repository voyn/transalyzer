function [bad_flag event_info event_more_info events_minmax voltage fsample timestep localtracetime peak_detection_factor sigma_from_detect_GUI readextrapoints totaleventsfound totaltimetrace rawunfilteredtrace_vect rawtrace_vect] = cmd_load_transalyzer_data(pathname, filename, verbose, check_load_rawtrace)

% this will load a Transalyzer standard dataset into the workspace

% example usage (use a tidle ~ if in place of unwanted output variables):
% filename = 'analysis_events';
% pathname = 'C:\Calin\Data\NanoPore\2013\2013-02-20\100mV\~transalyzed-2000Hz\';
% [bad_flag_low, event_info_low, event_more_info_low, ~, ~, ~, timestep, ~, ~, ~, ~, ~, ~, ~, rawtrace_vect_low] = load_transalyzer_data(pathname, filename, 1, 1);
%

% assign temporary values to output variables
bad_flag = 0;
event_info = -1;
event_more_info = -1;
events_minmax = -1;
voltage = -1;
fsample = -1;
timestep = -1;
localtracetime = -1;
peak_detection_factor = -1;
sigma_from_detect_GUI = -1;
readextrapoints = -1;
totaleventsfound = -1;
totaltimetrace = -1;
rawunfilteredtrace_vect = -1;
rawtrace_vect= -1;

if verbose == 1
    disp('Opening events file..')
end

fid = fopen(strcat(pathname,filesep,filename)); 

if fid == -1
    disp(strcat(filename,' could not be loaded!!')) % not good
    bad_flag = 1;
else
    % grab the event info
    event_info = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f')); 
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
    fclose(fid);
end

    
if verbose == 1 && bad_flag == 0
    disp('Opening events file...done...opening analysis_extra_info.txt...')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read analysis_extra_info.txt
% 
filenameextrainfo = strcat(pathname, filesep,'analysis_extra_info.txt');
if isequal(exist(filenameextrainfo),2) % if the extra info file exists
    % open it
    fid = fopen(filenameextrainfo); 
    if fid == -1
        disp(strcat(filenameextrainfo,' could not be loaded!!'))
        bad_flag = 1;
    else
        event_more_info = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f')); % get the extra info
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
        fclose(fid);
    end
end
if verbose == 1 && bad_flag == 0
    disp('Opening analysis_extra_info.txt...done...opening analysis_events_min_max...')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in analysis_events_min_max
filenameminmaxinfo = strcat(pathname, filesep,'analysis_events_min_max');
if isequal(exist(filenameminmaxinfo),2) % if the extra info file exists
    % open it
    fid = fopen(filenameminmaxinfo); 
    if fid == -1 % catch problems
        disp(strcat(filenameminmaxinfo,' could not be loaded!!'))
        bad_flag = 1;
    else
        % get the data
        events_minmax = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f %f')); 
        % this is for backwards compatibility
        events_minmax(:,8) = abs(events_minmax(:,8)); 
        % this is for backwards compatibility
        events_minmax(:,10) = abs(events_minmax(:,10)); 
        % select everything
        fclose(fid);
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

    
end
if verbose == 1 && bad_flag == 0
    disp('Opening analysis_events_min_max...done...opening analysis_parameters.txt...')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in analysis_parameters.txt
filenameparam = strcat(pathname, filesep,'analysis_parameters.txt');
if isequal(exist(filenameparam),2) % if the paramter file exists
    fid = fopen(filenameparam);
    if fid == -1 % catch problem
        disp(strcat(filenameparam,' could not be loaded!!'))
        bad_flag = 1;
    else % file loaded successfully
        line = fgetl(fid); % skip Transalyzer Analysis
        line = fgetl(fid); % skip version
        line = fgetl(fid); % skip Analysis Date
        line = fgetl(fid); % date, not used anywhere
        line = fgetl(fid); % skip
        voltage = str2double(fgetl(fid)); % Voltage (mV)
        line = fgetl(fid); % skip
        fsample = str2double(fgetl(fid)); % Sampling Frequency (Hz)
        timestep = 1/fsample;
        line = fgetl(fid); % skip
        localtracetime = str2double(fgetl(fid)); % Local Trace Time (s)
        line = fgetl(fid); % skip
        peak_detection_factor = str2double(fgetl(fid)); % Peak Det Factor
        line = fgetl(fid); % skip
        sigma_from_detect_GUI = str2double(fgetl(fid)); % Sigma (nA)
        line = fgetl(fid); % skip
        readextrapoints = str2double(fgetl(fid)); % Extra Points (points)
        line = fgetl(fid); % skip
        totaleventsfound = str2double(fgetl(fid)); % Total Events Found
        line = fgetl(fid); % skip
        totaltimetrace = str2double(fgetl(fid)); % Total Time
        fclose(fid); % close the file
    end
    
end
if verbose == 1 && bad_flag == 0
    disp('Opening analysis_parameters.txt...done...opening ...analysis_unfiltered_rawtrace')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in unfiltered current trace file (rawtrace)
if check_load_rawtrace == 1 % only load rawtrace if user wants to
    filenameunfilttrace = strcat(pathname, filesep,'analysis_unfiltered_rawtrace');
    if isequal(exist(filenameunfilttrace),2) % check file exists
        fid = fopen(filenameunfilttrace); % open it
        if fid == -1 % catch problems
            disp(strcat(filenameunfilttrace,' could not be loaded!!'))
            bad_flag = 1;
        else
            rawunfilteredtrace_vect = cell2mat(textscan(fid, '%f')); % read it in
            fclose(fid);
        end
        
    end
end
if verbose == 1 && bad_flag == 0
    disp('Opening analysis_unfiltered_rawtrace...done...opening ...rawtrace')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in current trace file (rawtrace)
if check_load_rawtrace == 1 % if user want to load rawtrace
    fid = fopen(strcat(pathname, filesep,'analysis_rawtrace'));
    if fid == -1 % catch problems
            disp(strcat(filename,' could not be loaded!!'))
            bad_flag = 1;
    else
        % read in rawtrace file
        rawtrace_vect = cell2mat(textscan(fid, '%f')); 
        fclose(fid);
    end
    
end
if verbose == 1 && bad_flag == 0
    disp('Opening rawtrace...done.')
end

