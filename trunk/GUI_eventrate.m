function varargout = GUI_eventrate(varargin)
% GUI_EVENTRATE M-file for GUI_eventrate.fig
%      GUI_EVENTRATE, by itself, creates a new GUI_EVENTRATE or raises the existing
%      singleton*.
%
%      H = GUI_EVENTRATE returns the handle to a new GUI_EVENTRATE or the handle to
%      the existing singleton*.
%
%      GUI_EVENTRATE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_EVENTRATE.M with the given input arguments.
%
%      GUI_EVENTRATE('Property','Value',...) creates a new GUI_EVENTRATE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_eventrate_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_eventrate_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_eventrate

% Last Modified by GUIDE v2.5 07-Mar-2012 19:21:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_eventrate_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_eventrate_OutputFcn, ...
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


% --- Executes just before GUI_eventrate is made visible.
function GUI_eventrate_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to GUI_eventrate (see VARARGIN)

    % Choose default command line output for GUI_eventrate
    handles.output = hObject;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isappdata(0,'HandleMainGUI')
        HandleMainGUI=getappdata(0,'HandleMainGUI');
        if isappdata(HandleMainGUI,'localtracetime') % if the variables were passed from the parent GUI

            % assign to local variables
            handles.pathname=getappdata(HandleMainGUI,'pathname');
            dwelltimes=getappdata(HandleMainGUI,'dwellvector');
            localtracetime=getappdata(HandleMainGUI,'localtracetime');
            handles.localtracetime = localtracetime;
            set(handles.text_localtime,'String',num2str(localtracetime))

            totaltimetrace=getappdata(HandleMainGUI,'totaltimetrace');
            set(handles.text_totaltime,'String',num2str(totaltimetrace))

            totaleventsfound=getappdata(HandleMainGUI,'totaleventsfound');
            set(handles.text_numcrude,'String',num2str(totaleventsfound))

            event_more_info_temp=getappdata(HandleMainGUI,'event_more_info_temp');
   
            
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

            % total number of event found, after sorting
            [numfine co] = size(event_more_info_temp);
            set(handles.text_numfine,'String',num2str(numfine))

            % event rate without sorting assuming uniform distribution
            er1 = totaleventsfound/totaltimetrace;
            set(handles.text_er1,'String',num2str(er1))

            % event rate with sorting assuming uniform distribution
            er2 = numfine/totaltimetrace;
            set(handles.text_er2,'String',num2str(er2))

            % number of files analyzed
            filenumbers = unique(event_more_info_temp(:,1), 'rows');
            %display(num2str(size(event_more_info_temp)))

            tracecounter = 0;
            averagecounter = 1;
            alleventcount=1;
            er3=0;
            er4=0;
            handles.er3_temp = 0;

            % go through each file or segment
            for filecounter=1:length(filenumbers)

                currentfile = filenumbers(filecounter);

                % get all of the indexes for event in this file
                eventindexes = find(event_more_info_temp(:,1)==currentfile);

                % number of events in this file
                handles.er3_temp(filecounter) = length(eventindexes);

                % get the number of the file
                handles.er3_filenumber(filecounter) = event_more_info_temp(eventindexes(1),1);

                if length(eventindexes) > 1 % if there is more than 1 event

                    clear averagetime

                    for eventcount = 1:length(eventindexes)-1 % for all events

                        % get differences between starting times of successive
                        % events
                        averagetime(eventcount)=event_more_info_temp(eventindexes(eventcount+1),3)-event_more_info_temp(eventindexes(eventcount),3);

                        % master vector of all event times
                        handles.timebetweeneventsall(alleventcount) = averagetime(eventcount);

                        % increment master event vector index
                        alleventcount = alleventcount + 1;
                    end
                    % get the mean for each segment
                    averagetimefile(averagecounter) = mean(averagetime);
                    averagecounter = averagecounter + 1;
                end

            end



%             for fillincounter = 1:max(filenumbers)
% 
%                 % get indexes of events for each file
%                 indexoffilinfile = find(handles.er3_filenumber==fillincounter);
% 
%                 if isempty(indexoffilinfile) % if a file had no events set ER as 0
%                     handles.number_of_events_in_each_segment(fillincounter) = 0;
%                 else % if a file had events set ER as the number of events
%                     handles.number_of_events_in_each_segment(fillincounter) = handles.er3_temp(indexoffilinfile);
%                 end
% 
%             end

            % divide by the time in each trace to get the rate
            handles.number_of_events_in_each_segment = handles.er3_temp./localtracetime;

            set(handles.edit_max_time,'String',num2str(max(handles.timebetweeneventsall)))

            er3 = mean(handles.number_of_events_in_each_segment); % average of all the rates
            er3_std = std(handles.number_of_events_in_each_segment); % error using std

            set(handles.text_er3,'String',num2str(er3))
            set(handles.text_er3_std,'String',num2str(er3_std))

            % Mean( Mean( local times ) )
            er4 = 1/(mean(averagetimefile));
            set(handles.text_er4,'String',num2str(er4))

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%% new stuff

            % here we calculate the time between events, on a global experimental time scale
            eventstarttimediff_experimental_temp = event_more_info_temp(:,4) - circshift(event_more_info_temp(:,4),[1 1]);

            eventstarttimediff_experimental_temp_include_dwell = event_more_info_temp(:,4) - circshift(event_more_info_temp(:,4),[1 1])-dwelltimes./1000;
            
            handles.eventstarttimediff_experimental = eventstarttimediff_experimental_temp(2:end); % discard first element
            eventstarttimediff_experimental_temp_include_dwell = eventstarttimediff_experimental_temp_include_dwell(2:end); % discard first element
            
            assignin('base', 'time_between_events_ms', eventstarttimediff_experimental_temp_include_dwell.*1000)
            assignin('base', 'events_dwell_ms', dwelltimes)
            
            % the largest time
            set(handles.edit_max_time2,'String',num2str(max(handles.eventstarttimediff_experimental)))

            er5 = 1/(mean(handles.eventstarttimediff_experimental));

            set(handles.text_er5,'String',num2str(er5))

            alphavalue = str2num(get(handles.edit_alpha_CI,'String'));

            % maximum time between events, local
            maxtime = str2num(get(handles.edit_max_time,'String'));
            % maximum time between events, global
            maxtime2 = str2num(get(handles.edit_max_time2,'String'));

            % resort using max time, local
            sortedtimebetweeneventsall_local = handles.timebetweeneventsall(find(handles.timebetweeneventsall <= maxtime));
            % resort using max time, all
            sortedtimebetweeneventsall_all = handles.eventstarttimediff_experimental(find(handles.eventstarttimediff_experimental <= maxtime2));

            % do the exponential fits

            % for local segments
            [meanval_local_time, confint_local_time] = expfit(sortedtimebetweeneventsall_local,alphavalue);
            set(handles.edit_expfitlocal_diffavg,'String',num2str(1/meanval_local_time))

            % for global
            [meanval_all_time, confint_all_time] = expfit(sortedtimebetweeneventsall_all,alphavalue);
            set(handles.edit_expfitall_diffavg,'String',num2str(1/meanval_all_time))

            % error using confidence interval
            set(handles.edit_expfit_error_high,'String',num2str(1/confint_all_time(1)-1/meanval_all_time))
            set(handles.edit_expfit_error_low,'String',num2str(1/meanval_all_time-1/confint_all_time(2)))

            plot_eventrate_vs_experimental_time(handles, 'old')

            plot_local_delta_time_and_fit(handles,'old',hObject)
            plot_all_delta_time_and_fit(handles,'old',hObject)

        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes GUI_eventrate wait for user response (see UIRESUME)
    % uiwait(handles.figure1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot a histogram of the times between successive events using times
% between files as well
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_all_delta_time_and_fit(handles, location, hObject)
    if strcmp(location, 'old')
        cla(handles.axeshist2,'reset')
        axes(handles.axeshist2)
    end
    if strcmp(location, 'new')
        %ylimits = get(handles.axeshist,'YLim');
        %xlimits = get(handles.axeshist,'XLim');
        figure()
    end
    %
    bins = str2num(get(handles.text_bins,'String'));
    
    % maximum time between events for plotting histogram
    maxtime2 = str2num(get(handles.edit_max_time2,'String'));
    
    % resort using max time
    handles.sorted_time_between_all_events_global = handles.eventstarttimediff_experimental(find(handles.eventstarttimediff_experimental <= maxtime2));
    %display(num2str(min(sortedtimebetweeneventsall2)))

    skip_first_X_bins = str2num(get(handles.edit_skip_first_bins,'String'));
    
    % do an exponential fit
    if (~get(handles.radio_lin_fit,'Value') & get(handles.radio_exp_fit,'Value'))
        
        % generate histogram
        [n_temp,xout_temp] = histnorm(handles.sorted_time_between_all_events_global,bins);
        
        hold on

        % skip first few bins in the fitting
        xout=xout_temp(skip_first_X_bins+1:end);
        n=n_temp(skip_first_X_bins+1:end);
        
        % calculate the location of the top edge
        xout_topedge = xout+(xout(2)-xout(1))*0.5.*ones(size(xout));
        
        histnorm(handles.sorted_time_between_all_events_global,bins)


        fit_var_type = fittype('exp1'); % a*exp(b*x)
        % Fit this model using new data
        fit_var_startingpoints = [1 -1];
        [fit_cf,fit_gof] = fit(xout',n',fit_var_type,'Startpoint',fit_var_startingpoints);
        % [fit_cf,fit_gof] = fit(xout',n',fit_var_type);
        
        plot(fit_cf,'r')
        legend off
        
        fit_mean_calc_time = -1/fit_cf.b;

        fit_mean_calc_rate = 1/fit_mean_calc_time;
        
        alphavalue = str2num(get(handles.edit_alpha_CI,'String'));
        ci = confint(fit_cf, 1-alphavalue); % 90% default
        %display(num2str(ci))
        confint_all_time = -1./[ci(1,2) ci(2,2)];

        
    elseif (get(handles.radio_lin_fit,'Value') & ~get(handles.radio_exp_fit,'Value')) % linear fit to log(n)
        
        [n_temp,xout_temp] = histnorm(handles.sorted_time_between_all_events_global,bins);
        
        % skip first few bins in the fitting
        xout = xout_temp(skip_first_X_bins+1:end);
        n = n_temp(skip_first_X_bins+1:end);
        
        n_log = log(n);
        n_log_clean = n_log(~isinf(n_log));
        xout_clean = xout(~isinf(n_log));
        
        plot(xout_clean,n_log_clean,'.b')
        hold on
        
        fit_var_type = fittype('poly1'); % f(x) = p1*x + p2
        % Fit this model using new data
        
        [fit_cf,fit_gof] = fit(xout_clean',n_log_clean',fit_var_type);
        
        
        plot(fit_cf,'r')
        legend off
        
        fit_mean_calc_time = -1/fit_cf.p1;

        fit_mean_calc_rate = 1/fit_mean_calc_time;
        
        alphavalue = str2num(get(handles.edit_alpha_CI,'String'));
        ci = confint(fit_cf, 1-alphavalue); % 90% default
        %display(num2str(ci))
        confint_all_time = -1./[ci(1,1) ci(2,1)];

        
    end
    
    confint_all_freq = 1./confint_all_time;
    diff_all = [confint_all_freq(1)-fit_mean_calc_rate fit_mean_calc_rate-confint_all_freq(2)];
    
    
    set(handles.edit_expfitall_mean,'String',num2str(fit_mean_calc_rate))

	set(handles.edit_expfitall_CIlow,'String',num2str(confint_all_freq(2)))
	set(handles.edit_expfitall_CIhigh,'String',num2str(confint_all_freq(1)))

    set(handles.edit_expfitall_difflow,'String',num2str(diff_all(2)))
	set(handles.edit_expfitall_diffhigh,'String',num2str(diff_all(1)))

    
    xlabel('Time (sec)')
    ylabel('Counts')
    if strcmp(location, 'old')
        %xlimits = get(handles.plot2,'XLim');
        %ylimits = get(handles.plot2,'YLim');
        %axis([0 (sampling_rate*(proteinlocation(1,1)+20)) 0 maxsig])

    end
    
    % Update handles structure
    guidata(hObject, handles);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the event rate as a function of experimental time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_eventrate_vs_experimental_time(handles, location)
    if strcmp(location, 'old')
        cla(handles.axeseventtime,'reset')
        axes(handles.axeseventtime)
    end
    if strcmp(location, 'new')
        ylimits = get(handles.axeseventtime,'YLim');
        xlimits = get(handles.axeseventtime,'XLim');
        figure()
    end
    
    timevector = handles.localtracetime:handles.localtracetime:handles.localtracetime*max(length(handles.number_of_events_in_each_segment));
    timevector = timevector - 0.5*handles.localtracetime.*ones(size(handles.number_of_events_in_each_segment));
    eventratevector = handles.number_of_events_in_each_segment./handles.localtracetime;
    plot(timevector, handles.number_of_events_in_each_segment, '-r')
    hold on
    plot(timevector, handles.number_of_events_in_each_segment, '.b','MarkerSize',15)
    %axis([0 (sampling_rate*(proteinlocation(1,1)+20)) 0 maxsig])
    xlabel('Time (s)')
    ylabel('Event Rate (Hz)')
    if strcmp(location, 'old')
        %xlimits = get(handles.plot2,'XLim');
        %ylimits = get(handles.plot2,'YLim');
        %axis([0 (sampling_rate*(proteinlocation(1,1)+20)) 0 maxsig])

    end
    if location == 'sav'
        dir_str = strcat(handles.pathname,'plots', filesep); % in the plots dir
        if(~exist(dir_str, 'dir'))
            mkdir(dir_str);
        end
        file_str = strcat(dir_str,'plot_data_eventrate_vs_experimental_time.dat'); % save dat file
        dataoutput = [timevector, handles.number_of_events_in_each_segment];
        dlmwrite(file_str, dataoutput, 'delimiter', '\t', 'precision', '%.8f','newline', 'pc');
    end
    if location == 'new'
        axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
    end    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot a histogram of the times between successive events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_local_delta_time_and_fit(handles, location, hObject)
    if strcmp(location, 'old')
        cla(handles.axeshist,'reset')
        axes(handles.axeshist)
    end
    if strcmp(location, 'new')
        %ylimits = get(handles.axeshist,'YLim');
        %xlimits = get(handles.axeshist,'XLim');
        figure()
    end
    bins = str2num(get(handles.text_bins,'String'));
    
    % maximum time between events for plotting histogram
    maxtime = str2num(get(handles.edit_max_time,'String'));
    
    % resort using max time
    handles.sortedtimebetweeneventsall_local = handles.timebetweeneventsall(find(handles.timebetweeneventsall <= maxtime));
        
    skip_first_X_bins = str2num(get(handles.edit_skip_first_bins,'String'));
    
    if (~get(handles.radio_lin_fit,'Value') & get(handles.radio_exp_fit,'Value')) %exp fit
    
        % generate histogram
        [n_temp,xout_temp] = histnorm(handles.sortedtimebetweeneventsall_local,bins);
        
        histnorm(handles.sortedtimebetweeneventsall_local,bins)
        
        hold on
              
        % skip first few bins in the fitting
        xout=xout_temp(skip_first_X_bins+1:end);
        n=n_temp(skip_first_X_bins+1:end);

        fit_var_type = fittype('exp1'); % a*exp(b*x)
        
        % Fit this model using new data
        fit_var_startingpoints = [1 -1];
        [fit_cf,fit_gof] = fit(xout',n',fit_var_type,'Startpoint',fit_var_startingpoints);
        
        plot(fit_cf,'r')
        legend off
        
        fit_mean_calc_time = -1/fit_cf.b;

        fit_mean_calc_rate = 1/fit_mean_calc_time;

        alphavalue = str2num(get(handles.edit_alpha_CI,'String'));
        ci = confint(fit_cf, 1-alphavalue); % 90% default
        %display(num2str(ci))
        confint_all_time = -1./[ci(1,2) ci(2,2)];
                
        
    elseif (get(handles.radio_lin_fit,'Value') & ~get(handles.radio_exp_fit,'Value')) % linear fit to log(n)
        
        % generate histogram
        [n_temp,xout_temp] = histnorm(handles.sortedtimebetweeneventsall_local,bins);
                     
        % skip first few bins in the fitting
        xout=xout_temp(skip_first_X_bins+1:end);
        n=n_temp(skip_first_X_bins+1:end);
        
        n_log = log(n);
        n_log_clean = n_log(~isinf(n_log));
        xout_clean = xout(~isinf(n_log));
        
        plot(xout_clean,n_log_clean,'.b')
        hold on

        fit_var_type = fittype('poly1'); % f(x) = p1*x + p2
        
        % Fit this model using new data
        fit_var_startingpoints = [1 -1];
        [fit_cf,fit_gof] = fit(xout_clean',n_log_clean',fit_var_type);
        
        plot(fit_cf,'r')
        legend off
        
        fit_mean_calc_time = -1/fit_cf.p1;

        fit_mean_calc_rate = 1/fit_mean_calc_time;

        alphavalue = str2num(get(handles.edit_alpha_CI,'String'));
        ci = confint(fit_cf, 1-alphavalue); % 90% default
        %display(num2str(ci))
        confint_all_time = -1./[ci(1,1) ci(2,1)];
        
            
            
        
    end
    
    confint_all_freq = 1./confint_all_time;
    diff_all = [confint_all_freq(1)-fit_mean_calc_rate fit_mean_calc_rate-confint_all_freq(2)];

    set(handles.edit_expfitlocal_mean,'String',num2str(fit_mean_calc_rate))

	set(handles.edit_expfitlocal_CIlow,'String',num2str(confint_all_freq(2)))
	set(handles.edit_expfitlocal_CIhigh,'String',num2str(confint_all_freq(1)))

    set(handles.edit_expfitlocal_difflow,'String',num2str(diff_all(2)))
	set(handles.edit_expfitlocal_diffhigh,'String',num2str(diff_all(1)))
	%set(handles.edit_expfitlocal_diffavg,'String',num2str(mean(diff_all)))
    
    xlabel('Time (sec)')
    ylabel('Counts')
    if strcmp(location, 'old')
        %xlimits = get(handles.plot2,'XLim');
        %ylimits = get(handles.plot2,'YLim');
        %axis([0 (sampling_rate*(proteinlocation(1,1)+20)) 0 maxsig])

    end
    % if location == 'new'
    %     axis([xlimits(1) xlimits(2) ylimits(1) ylimits(2)]);
    % end
    
    % Update handles structure
    guidata(hObject, handles);
    
% --- Outputs from this function are returned to the command line.
function varargout = GUI_eventrate_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_exit.
function button_exit_Callback(hObject, eventdata, handles)
% hObject    handle to button_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(gcf)


function text_localtime_Callback(hObject, eventdata, handles)
% hObject    handle to text_localtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_localtime as text
%        str2double(get(hObject,'String')) returns contents of text_localtime as a double


% --- Executes during object creation, after setting all properties.
function text_localtime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_localtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_totaltime_Callback(hObject, eventdata, handles)
% hObject    handle to text_totaltime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_totaltime as text
%        str2double(get(hObject,'String')) returns contents of text_totaltime as a double


% --- Executes during object creation, after setting all properties.
function text_totaltime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_totaltime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_numcrude_Callback(hObject, eventdata, handles)
% hObject    handle to text_numcrude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_numcrude as text
%        str2double(get(hObject,'String')) returns contents of text_numcrude as a double


% --- Executes during object creation, after setting all properties.
function text_numcrude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_numcrude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_numfine_Callback(hObject, eventdata, handles)
% hObject    handle to text_numfine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_numfine as text
%        str2double(get(hObject,'String')) returns contents of text_numfine as a double


% --- Executes during object creation, after setting all properties.
function text_numfine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_numfine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_er1_Callback(hObject, eventdata, handles)
% hObject    handle to text_er1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_er1 as text
%        str2double(get(hObject,'String')) returns contents of text_er1 as a double


% --- Executes during object creation, after setting all properties.
function text_er1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_er1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_er2_Callback(hObject, eventdata, handles)
% hObject    handle to text_er2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_er2 as text
%        str2double(get(hObject,'String')) returns contents of text_er2 as a double


% --- Executes during object creation, after setting all properties.
function text_er2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_er2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_er3_Callback(hObject, eventdata, handles)
% hObject    handle to text_er3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_er3 as text
%        str2double(get(hObject,'String')) returns contents of text_er3 as a double


% --- Executes during object creation, after setting all properties.
function text_er3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_er3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_er4_Callback(hObject, eventdata, handles)
% hObject    handle to text_er4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_er4 as text
%        str2double(get(hObject,'String')) returns contents of text_er4 as a double


% --- Executes during object creation, after setting all properties.
function text_er4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_er4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in button_openplot.
function button_openplot_Callback(hObject, eventdata, handles)
% hObject    handle to button_openplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_local_delta_time_and_fit(handles,'new',hObject)
% get updated handles
handles=guidata(gcbo);
plot_eventrate_vs_experimental_time(handles, 'new')
plot_all_delta_time_and_fit(handles,'new',hObject)
guidata(gcbo, handles);


% --- Executes on button press in button_replot.
function button_replot_Callback(hObject, eventdata, handles)
% hObject    handle to button_replot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_local_delta_time_and_fit(handles,'old',hObject)
% get updated handles
handles=guidata(gcbo);
plot_eventrate_vs_experimental_time(handles, 'old')
plot_all_delta_time_and_fit(handles,'old',hObject)
guidata(gcbo, handles);

function text_bins_Callback(hObject, eventdata, handles)
% hObject    handle to text_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_bins as text
%        str2double(get(hObject,'String')) returns contents of text_bins as a double


% --- Executes during object creation, after setting all properties.
function text_bins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in button_dump.
function button_dump_Callback(hObject, eventdata, handles)
    % hObject    handle to button_dump (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % TODO probably stuff missing here
    
    if isfield(handles,'event_info')
        assignin('base', 'event_info', handles.event_info)
    end
    if isfield(handles,'localtracetime') && isfield(handles,'number_of_events_in_each_segment')
        timevector = handles.localtracetime:handles.localtracetime:handles.localtracetime*max(length(handles.number_of_events_in_each_segment));
        timevector = timevector - 0.5*handles.localtracetime.*ones(size(handles.number_of_events_in_each_segment));
        assignin('base', 'timevector', timevector')
        assignin('base', 'number_of_events_in_each_segment', handles.number_of_events_in_each_segment')
    end
    if isfield(handles,'event_info')
        assignin('base', 'event_info', handles.event_info)
    end
    if isfield(handles,'sortedtimebetweeneventsall_local')
        assignin('base', 'sortedtimebetweeneventsall_local', handles.sortedtimebetweeneventsall_local)
    end
    if isfield(handles,'sorted_time_between_all_events_global')
        assignin('base', 'sorted_time_between_all_events_global', handles.sorted_time_between_all_events_global)
    end
    



% --- Executes on button press in button_saveplots.
function button_saveplots_Callback(hObject, eventdata, handles)
    % hObject    handle to button_saveplots (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    dir_str = strcat(handles.pathname,'plots', filesep);
    if(~exist(dir_str, 'dir'))
        mkdir(dir_str);
    end

    newFig = figure;
    plot_local_delta_time_and_fit(handles, 'sav')
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_event_separation_hist.fig'), 'fig');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_event_separation_hist_large.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_event_separation_hist_small.png'), 'png');
    close(newFig)

    newFig = figure;
    plot_eventrate_vs_experimental_time(handles, 'sav')
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_eventrate.fig'), 'fig');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_eventrate.png'), 'png');
    set(newFig, 'PaperPositionMode', 'auto');
    saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_eventrate.png'), 'png');
    close(newFig)

    %TODO add the other plots





function edit_max_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max_time as text
%        str2double(get(hObject,'String')) returns contents of edit_max_time as a double


% --- Executes during object creation, after setting all properties.
function edit_max_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function text_er5_Callback(hObject, eventdata, handles)
% hObject    handle to text_er5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_er5 as text
%        str2double(get(hObject,'String')) returns contents of text_er5 as a double


% --- Executes during object creation, after setting all properties.
function text_er5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_er5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_er5_std_Callback(hObject, eventdata, handles)
% hObject    handle to text_er5_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_er5_std as text
%        str2double(get(hObject,'String')) returns contents of text_er5_std as a double


% --- Executes during object creation, after setting all properties.
function text_er5_std_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_er5_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_er1_std_Callback(hObject, eventdata, handles)
% hObject    handle to text_er1_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_er1_std as text
%        str2double(get(hObject,'String')) returns contents of text_er1_std as a double


% --- Executes during object creation, after setting all properties.
function text_er1_std_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_er1_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_er2_std_Callback(hObject, eventdata, handles)
% hObject    handle to text_er2_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_er2_std as text
%        str2double(get(hObject,'String')) returns contents of text_er2_std as a double


% --- Executes during object creation, after setting all properties.
function text_er2_std_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_er2_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_er3_std_Callback(hObject, eventdata, handles)
% hObject    handle to text_er3_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_er3_std as text
%        str2double(get(hObject,'String')) returns contents of text_er3_std as a double


% --- Executes during object creation, after setting all properties.
function text_er3_std_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_er3_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_er4_std_Callback(hObject, eventdata, handles)
% hObject    handle to text_er4_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_er4_std as text
%        str2double(get(hObject,'String')) returns contents of text_er4_std as a double


% --- Executes during object creation, after setting all properties.
function text_er4_std_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_er4_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_openplot2.
function button_openplot2_Callback(hObject, eventdata, handles)
% hObject    handle to button_openplot2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_replot2.
function button_replot2_Callback(hObject, eventdata, handles)
% hObject    handle to button_replot2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_max_time2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max_time2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max_time2 as text
%        str2double(get(hObject,'String')) returns contents of edit_max_time2 as a double


% --- Executes during object creation, after setting all properties.
function edit_max_time2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max_time2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_expfitall_mean_Callback(hObject, eventdata, handles)
% hObject    handle to edit_expfitall_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_expfitall_mean as text
%        str2double(get(hObject,'String')) returns contents of edit_expfitall_mean as a double


% --- Executes during object creation, after setting all properties.
function edit_expfitall_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_expfitall_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_expfitall_CIlow_Callback(hObject, eventdata, handles)
% hObject    handle to edit_expfitall_CIlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_expfitall_CIlow as text
%        str2double(get(hObject,'String')) returns contents of edit_expfitall_CIlow as a double


% --- Executes during object creation, after setting all properties.
function edit_expfitall_CIlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_expfitall_CIlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_expfitall_CIhigh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_expfitall_CIhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_expfitall_CIhigh as text
%        str2double(get(hObject,'String')) returns contents of edit_expfitall_CIhigh as a double


% --- Executes during object creation, after setting all properties.
function edit_expfitall_CIhigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_expfitall_CIhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_alpha_CI_Callback(hObject, eventdata, handles)
% hObject    handle to edit_alpha_CI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_alpha_CI as text
%        str2double(get(hObject,'String')) returns contents of edit_alpha_CI as a double


% --- Executes during object creation, after setting all properties.
function edit_alpha_CI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_alpha_CI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_expfitlocal_mean_Callback(hObject, eventdata, handles)
% hObject    handle to edit_expfitlocal_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_expfitlocal_mean as text
%        str2double(get(hObject,'String')) returns contents of edit_expfitlocal_mean as a double


% --- Executes during object creation, after setting all properties.
function edit_expfitlocal_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_expfitlocal_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_expfitlocal_CIlow_Callback(hObject, eventdata, handles)
% hObject    handle to edit_expfitlocal_CIlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_expfitlocal_CIlow as text
%        str2double(get(hObject,'String')) returns contents of edit_expfitlocal_CIlow as a double


% --- Executes during object creation, after setting all properties.
function edit_expfitlocal_CIlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_expfitlocal_CIlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_expfitlocal_CIhigh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_expfitlocal_CIhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_expfitlocal_CIhigh as text
%        str2double(get(hObject,'String')) returns contents of edit_expfitlocal_CIhigh as a double


% --- Executes during object creation, after setting all properties.
function edit_expfitlocal_CIhigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_expfitlocal_CIhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_expfitall_difflow_Callback(hObject, eventdata, handles)
% hObject    handle to edit_expfitall_difflow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_expfitall_difflow as text
%        str2double(get(hObject,'String')) returns contents of edit_expfitall_difflow as a double


% --- Executes during object creation, after setting all properties.
function edit_expfitall_difflow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_expfitall_difflow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_expfitall_diffhigh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_expfitall_diffhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_expfitall_diffhigh as text
%        str2double(get(hObject,'String')) returns contents of edit_expfitall_diffhigh as a double


% --- Executes during object creation, after setting all properties.
function edit_expfitall_diffhigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_expfitall_diffhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_expfitlocal_difflow_Callback(hObject, eventdata, handles)
% hObject    handle to edit_expfitlocal_difflow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_expfitlocal_difflow as text
%        str2double(get(hObject,'String')) returns contents of edit_expfitlocal_difflow as a double


% --- Executes during object creation, after setting all properties.
function edit_expfitlocal_difflow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_expfitlocal_difflow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_expfitlocal_diffhigh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_expfitlocal_diffhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_expfitlocal_diffhigh as text
%        str2double(get(hObject,'String')) returns contents of edit_expfitlocal_diffhigh as a double


% --- Executes during object creation, after setting all properties.
function edit_expfitlocal_diffhigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_expfitlocal_diffhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_expfitall_diffavg_Callback(hObject, eventdata, handles)
% hObject    handle to edit_expfitall_diffavg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_expfitall_diffavg as text
%        str2double(get(hObject,'String')) returns contents of edit_expfitall_diffavg as a double


% --- Executes during object creation, after setting all properties.
function edit_expfitall_diffavg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_expfitall_diffavg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_expfitlocal_diffavg_Callback(hObject, eventdata, handles)
% hObject    handle to edit_expfitlocal_diffavg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_expfitlocal_diffavg as text
%        str2double(get(hObject,'String')) returns contents of edit_expfitlocal_diffavg as a double


% --- Executes during object creation, after setting all properties.
function edit_expfitlocal_diffavg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_expfitlocal_diffavg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_skip_first_bins_Callback(hObject, eventdata, handles)
% hObject    handle to edit_skip_first_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_skip_first_bins as text
%        str2double(get(hObject,'String')) returns contents of edit_skip_first_bins as a double


% --- Executes during object creation, after setting all properties.
function edit_skip_first_bins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_skip_first_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_expfit_error_high_Callback(hObject, eventdata, handles)
% hObject    handle to edit_expfit_error_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_expfit_error_high as text
%        str2double(get(hObject,'String')) returns contents of edit_expfit_error_high as a double


% --- Executes during object creation, after setting all properties.
function edit_expfit_error_high_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_expfit_error_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_expfit_error_low_Callback(hObject, eventdata, handles)
% hObject    handle to edit_expfit_error_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_expfit_error_low as text
%        str2double(get(hObject,'String')) returns contents of edit_expfit_error_low as a double


% --- Executes during object creation, after setting all properties.
function edit_expfit_error_low_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_expfit_error_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
