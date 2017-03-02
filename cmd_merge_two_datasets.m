% merge two analysis datasets
% full paths to directories, must end in file separator / or \
analysis_path_1 = '/Users/calin/Documents/phd/analysis_1/';
analysis_path_2 = '/Users/calin/Documents/phd/analysis_2/';
merged_set_path = '/Users/calin/Documents/phd/merged/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open all files
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read analysis_events 1
% 
disp('Opening events file 1 ...')
% open events file 1
fullfilepath = strcat(analysis_path_1,'analysis_events');
fid = fopen(fullfilepath); 
if fid == -1
    disp(strcat(fullfilepath,' could not be loaded!!')) % not good
else
    % grab the event info
    event_info_1 = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f')); 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read analysis_events 2
% 
disp('Opening events file 2 ...')
% open events file 1
fullfilepath = strcat(analysis_path_2,'analysis_events');
fid = fopen(fullfilepath); 
if fid == -1
    disp(strcat(fullfilepath,' could not be loaded!!')) % not good
else
    % grab the event info
    event_info_2 = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f')); 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read analysis_extra_info.txt 1
% 
filenameextrainfo = strcat(analysis_path_1,'analysis_extra_info.txt');
if isequal(exist(filenameextrainfo),2) % if the extra info file exists
    % open it
    fid = fopen(filenameextrainfo); 
    if fid == -1
        disp(strcat(filenameextrainfo,' could not be loaded!!'))
    else
        event_more_info_1 = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f')); % get the extra info
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
% read analysis_extra_info.txt 2
% 
filenameextrainfo = strcat(analysis_path_2,'analysis_extra_info.txt');
if isequal(exist(filenameextrainfo),2) % if the extra info file exists
    % open it
    fid = fopen(filenameextrainfo); 
    if fid == -1
        disp(strcat(filenameextrainfo,' could not be loaded!!'))
    else
        event_more_info_2 = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f')); % get the extra info
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
% read in analysis_events_min_max 1
filenameminmaxinfo = strcat(analysis_path_1,'analysis_events_min_max');
if isequal(exist(filenameminmaxinfo),2) % if the extra info file exists
    % open it
    fid = fopen(filenameminmaxinfo); 
    if fid == -1 % catch problems
        disp(strcat(filenameminmaxinfo,' could not be loaded!!'))
    else
        % get the data
        events_minmax_1 = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f %f')); 
        % this is for backwards compatibility
        events_minmax_1(:,8) = abs(events_minmax_1(:,8)); 
        % this is for backwards compatibility
        events_minmax_1(:,10) = abs(events_minmax_1(:,10)); 
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
% read in analysis_events_min_max 2
filenameminmaxinfo = strcat(analysis_path_2,'analysis_events_min_max');
if isequal(exist(filenameminmaxinfo),2) % if the extra info file exists
    % open it
    fid = fopen(filenameminmaxinfo); 
    if fid == -1 % catch problems
        disp(strcat(filenameminmaxinfo,' could not be loaded!!'))
    else
        % get the data
        events_minmax_2 = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f %f')); 
        % this is for backwards compatibility
        events_minmax_2(:,8) = abs(events_minmax_2(:,8)); 
        % this is for backwards compatibility
        events_minmax_2(:,10) = abs(events_minmax_2(:,10)); 
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
% read in unfiltered current trace file 1
filenameunfilttrace = strcat(analysis_path_1,'analysis_unfiltered_rawtrace');
if isequal(exist(filenameunfilttrace),2) % check file exists
    fid = fopen(filenameunfilttrace); % open it
    if fid == -1 % catch problems
        disp(strcat(filenameunfilttrace,' could not be loaded!!'))
    else
        rawunfilteredtrace_vect_1 = cell2mat(textscan(fid, '%f')); % read it in
    end
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in unfiltered current trace file 2
filenameunfilttrace = strcat(analysis_path_2,'analysis_unfiltered_rawtrace');
if isequal(exist(filenameunfilttrace),2) % check file exists
    fid = fopen(filenameunfilttrace); % open it
    if fid == -1 % catch problems
        disp(strcat(filenameunfilttrace,' could not be loaded!!'))
    else
        rawunfilteredtrace_vect_2 = cell2mat(textscan(fid, '%f')); % read it in
    end
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in trace file (rawtrace) 1
filenametrace = strcat(analysis_path_1,'analysis_rawtrace');
fid = fopen(filenametrace);
if fid == -1 % catch problems
        disp(strcat(filenametrace,' could not be loaded!!'))
else
    % read in rawtrace file
    rawtrace_vect_1 = cell2mat(textscan(fid, '%f')); 
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in trace file (rawtrace) 2
filenametrace = strcat(analysis_path_2,'analysis_rawtrace');
fid = fopen(filenametrace);
if fid == -1 % catch problems
        disp(strcat(filenametrace,' could not be loaded!!'))
else
    % read in rawtrace file
    rawtrace_vect_2 = cell2mat(textscan(fid, '%f')); 
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge all variables
%

if (isequal(exist('rawtrace_vect_1'),1) && isequal(exist('rawtrace_vect_2'),1)) % check file exists
    % update the offset for merged rawtrace
    rawtrace_offset = length(rawtrace_vect_1);
    event_info_2(:,8) = event_info_2(:,8) + rawtrace_offset;
    event_info_2(:,9) = event_info_2(:,9) + rawtrace_offset;
    rawtrace_vect = vertcat(rawtrace_vect_1,rawtrace_vect_2);
end

event_info = vertcat(event_info_1,event_info_2);

if (isequal(exist('event_more_info_1'),1) && isequal(exist('event_more_info_2'),1)) % check file exists
    event_more_info = vertcat(event_more_info_1,event_more_info_2);
end

if (isequal(exist('events_minmax_1'),1) && isequal(exist('events_minmax_2'),1)) % check file exists
    events_minmax = vertcat(events_minmax_1,events_minmax_2);
end

if (isequal(exist('rawunfilteredtrace_vect_1'),1) && isequal(exist('rawunfilteredtrace_vect_2'),1)) % check file exists
    rawunfilteredtrace_vect = vertcat(rawunfilteredtrace_vect_1,rawunfilteredtrace_vect_2);
end

disp(strcat(num2str(size(event_info_1,1)), ' events in 1st set'))
disp(strcat(num2str(size(event_info_2,1)), ' events in 2nd set'))
disp(strcat(num2str(size(event_info,1)), ' events in merged set'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save merged variables
%

% write events file
dlmwrite(strcat(merged_set_path,'analysis_events'), event_info, 'delimiter', '\t', 'precision', '%.8f');

if isequal(exist('events_minmax'),1) % write minmax file
    dlmwrite(strcat(merged_set_path,'analysis_events_min_max'), events_minmax, 'delimiter', '\t', 'precision', '%.8f'); 
end

% write extra info file
dlmwrite(strcat(merged_set_path,'analysis_extra_info.txt'), event_more_info, 'delimiter', '\t', 'precision', '%.8f'); 

if isequal(exist('rawtrace_vect'),1)
    % write rawtrace file
    dlmwrite(strcat(merged_set_path,'analysis_rawtrace'), rawtrace_vect, 'delimiter', '\n', 'precision', '%.8f'); 
    if isequal(exist('rawunfilteredtrace_vect'),1)
        % write unfiltered rawtrace
        dlmwrite(strcat(merged_set_path,'analysis_unfiltered_rawtrace'), rawunfilteredtrace_vect, 'delimiter', '\n', 'precision', '%.8f'); 
    end
end

% copy ~detect.par
copyfile(strcat(analysis_path_1,'~detect.par'),strcat(merged_set_path,'~detect.par'))
% copy legacy file
copyfile(strcat(analysis_path_1,'analysis_parameters.txt'),strcat(merged_set_path,'analysis_parameters.txt'))

disp('Done writing out.')
