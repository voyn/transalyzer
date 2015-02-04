function varargout = GUI_noise(varargin)
% GUI_NOISE M-file for GUI_noise.fig
%      GUI_NOISE, by itself, creates a new GUI_NOISE or raises the existing
%      singleton*.
%
%      H = GUI_NOISE returns the handle to a new GUI_NOISE or the handle to
%      the existing singleton*.
%
%      GUI_NOISE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_NOISE.M with the given input arguments.
%
%      GUI_NOISE('Property','Value',...) creates a new GUI_NOISE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_noise_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_noise_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help GUI_noise

% Last Modified by GUIDE v2.5 15-Mar-2012 15:18:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_noise_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_noise_OutputFcn, ...
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


% --- Executes just before GUI_noise is made visible.
function GUI_noise_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_noise (see VARARGIN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flagin = 0;
% if isappdata(0,'HandleMainGUI')
%     HandleMainGUI=getappdata(0,'HandleMainGUI');
%     if isappdata(HandleMainGUI,'noisetrace')
%     	handles.filetype_in = getappdata(HandleMainGUI,'filetype');
%         handles.t1traceraw = getappdata(HandleMainGUI,'noisetrace');
%         handles.t1timestep = getappdata(HandleMainGUI,'noisetimestep');
%         flagin = 1;
%     end
% else
%     
% end
pwd;
currentFolder = pwd;
addpath(currentFolder)

% Choose default command line output for GUI_noise
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if flagin == 1
    set(handles.t1time,'String',strcat(num2str(handles.t1timestep*length(handles.t1traceraw)),' seconds long'))
    set(handles.statustext,'String','Trace 1 Passed in From GUI_detect')
    fun_t1plot(handles, 1) % trace came in
end

% UIWAIT makes GUI_noise wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_noise_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in buttonexit.
function buttonexit_Callback(hObject, eventdata, handles)
% hObject    handle to buttonexit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(gcf)



% --- Executes on button press in t1mult.
function t1mult_Callback(hObject, eventdata, handles)
% hObject    handle to t1mult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of t1mult


% --- Executes on button press in t1smooth.
function t1smooth_Callback(hObject, eventdata, handles)
% hObject    handle to t1smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of t1smooth


% --- Executes on button press in t1dc.
function t1dc_Callback(hObject, eventdata, handles)
% hObject    handle to t1dc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of t1dc


% --- Executes on button press in t1load.
function t1load_Callback(hObject, eventdata, handles)
    % hObject    handle to t1load (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    [handles.t1filename, handles.t1pathname] = uigetfile( ...
    {  '*',  'All Files (*.*)'}, ...
       'Pick a file');
    set(handles.t1fname,'String',handles.t1filename)
    switch get(handles.t1pop,'Value')   
        case 1 % Record Trace Bin
            [t1trace, t1time_vector, t1timestep, t1code] = readlabviewbinaries(strcat(handles.t1pathname,handles.t1filename));
            if get(handles.t1mult,'Value') == 1
                t1trace = handles.t1traceraw.*10;
            end
            t1trace = t1trace.*1000;
        case 2 % ATF file
            fid = fopen(strcat(handles.t1pathname,handles.t1filename));
            if fid > 0 % opened
                for i = 1:10 % header to be skipped
                    line = fgetl(fid);
                end
                fffnoise = fscanf(fid, '%f\t%f', [2 inf]);
                
                t1timestep = fffnoise(1,3)-fffnoise(1,2);
                t1trace = fffnoise(2,:); % in pA
            end
            fclose(fid);
        case 3 % Text file with two columns, same as above but without header
            fid = fopen(strcat(handles.t1pathname,handles.t1filename));
            if fid > 0
                fffnoise = fscanf(fid, '%f\t%f', [2 inf]);
                t1timestep = fffnoise(1,3)-fffnoise(1,2); %% column 1 is time
                t1trace = fffnoise(2,:).*1000; %% column 2 is current
            end
            fclose(fid);
        case 4 % ABF file
            segment_size_in_sec = str2num(get(handles.text_t1_duration,'String'));
            start_time = str2num(get(handles.text_t1_start,'String'));
                        
            [t1trace,si,lActualAcqLength]=abfload(strcat(handles.t1pathname,handles.t1filename),'start',start_time,'stop',segment_size_in_sec);
            
            t1timestep = si*10^(-6); % convert to sec
            t1trace = transpose(t1trace); % keep in pA, and switch to correct form
        otherwise
    end
    handles.t1timestep = t1timestep;
    handles.t1traceraw = t1trace;
    set(handles.t1time,'String',strcat(num2str(handles.t1timestep*length(handles.t1traceraw)),' seconds long'))
    set(handles.statustext,'String','Trace 1 File Loaded')
    %t1refresh(handles);
    guidata(gcbo, handles); % update handles


    
function t1fname_Callback(hObject, eventdata, handles)
% hObject    handle to t1fname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t1fname as text
%        str2double(get(hObject,'String')) returns contents of t1fname as a double


% --- Executes during object creation, after setting all properties.
function t1fname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t1fname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in t1plot.
function t1plot_Callback(hObject, eventdata, handles)
    % hObject    handle to t1plot (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    fun_t1plot(handles, 0)
    
function fun_t1plot(handles, flagin)
    t1trace = handles.t1traceraw;    
    
    if get(handles.t1dc,'Value') == 1
        t1trace = t1trace-mean(t1trace);
    end
    FsG = 1/(handles.t1timestep);
    [fGn PxxGn] = calcPSD(FsG, t1trace);
    if get(handles.t1smooth,'Value') == 1
        [smoothed_fGn smoothed_PxxGn] = smoothing(fGn, PxxGn);
        fGn = smoothed_fGn; 
        PxxGn = smoothed_PxxGn;
    end
    handles.t1fGn = fGn;
    handles.t1PxxGn = PxxGn;
    t1plotty(handles, 'here')
    
    %set(handles.statustext,'String',strcat('Pars_time = ', num2str(Pars_time),' and Pars_freq = ',num2str(Pars_freq),' and Pt/Pf = ',num2str(Ratio_Time_Freq)))
    guidata(gcbo, handles); % update handles
    
function t1plotty(handles, location)
    if get(handles.lockaxes,'Value') == 1
        lxmin = str2num(get(handles.xmin,'String'));
        lxmax = str2num(get(handles.xmax,'String'));
        lymin = str2num(get(handles.ymin,'String'));
        lymax = str2num(get(handles.ymax,'String'));
    end
    if strcmp(location, 'here')
        cla(handles.axes1,'reset')
        axes(handles.axes1)
    end
    if strcmp(location, 'new')
        figure()
    end
    loglog(handles.t1fGn,handles.t1PxxGn,'-r')
    xlabel('Frequency (Hz)')
    ylabel('PSD (nA^2/Hz)')  
    if get(handles.lockaxes,'Value') == 1
        axis([lxmin lxmax lymin lymax]);
    end
    xlimits = get(handles.axes1,'XLim');
    ylimits = get(handles.axes1,'YLim');
    
    set(handles.xmin,'String',num2str(xlimits(1)))
    set(handles.xmax,'String',num2str(xlimits(2)))
    set(handles.ymin,'String',num2str(ylimits(1)))
    set(handles.ymax,'String',num2str(ylimits(2)))



function xmin_Callback(hObject, eventdata, handles)
% hObject    handle to xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xmin as text
%        str2double(get(hObject,'String')) returns contents of xmin as a double
set(handles.statustext,'String','Setting Plot Xmin')
xlimits = get(handles.axes1,'XLim');
minvalue = str2double(get(handles.xmin,'String'));
set(handles.axes1,'XLim',[minvalue xlimits(2)]);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function xmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xmax_Callback(hObject, eventdata, handles)
% hObject    handle to xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xmax as text
%        str2double(get(hObject,'String')) returns contents of xmax as a double
set(handles.statustext,'String','Setting Plot Xmax')
xlimits = get(handles.axes1,'XLim');
maxvalue = str2double(get(handles.xmax,'String'));
set(handles.axes1,'XLim',[xlimits(1) maxvalue]);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function xmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ymin_Callback(hObject, eventdata, handles)
% hObject    handle to ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ymin as text
%        str2double(get(hObject,'String')) returns contents of ymin as a double
set(handles.statustext,'String','Setting Plot Ymin')
ylimits = get(handles.axes1,'YLim');
minvalue = str2double(get(handles.ymin,'String'));
set(handles.axes1,'YLim',[minvalue ylimits(2)]);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ymin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ymax_Callback(hObject, eventdata, handles)
% hObject    handle to ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ymax as text
%        str2double(get(hObject,'String')) returns contents of ymax as a double
set(handles.statustext,'String','Setting Plot Ymax')
ylimits = get(handles.axes1,'YLim');
maxvalue = str2double(get(handles.ymax,'String'));
set(handles.axes1,'YLim',[ylimits(1) maxvalue]);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ymax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in t1enable.
function t1enable_Callback(hObject, eventdata, handles)
% hObject    handle to t1enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of t1enable


% --- Executes on selection change in t1pop.
function t1pop_Callback(hObject, eventdata, handles)
% hObject    handle to t1pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns t1pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from t1pop


% --- Executes during object creation, after setting all properties.
function t1pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t1pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in t2mult.
function t2mult_Callback(hObject, eventdata, handles)
% hObject    handle to t2mult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of t2mult


% --- Executes on button press in t2dc.
function t2dc_Callback(hObject, eventdata, handles)
% hObject    handle to t2dc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of t2dc


% --- Executes on button press in t2load.
function t2load_Callback(hObject, eventdata, handles)
% hObject    handle to t2load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [handles.t2filename, handles.t2pathname] = uigetfile( ...
    {  '*',  'All Files (*.*)'}, ...
       'Pick a file');
    set(handles.t2fname,'String',handles.t2filename)
    switch get(handles.t2pop,'Value')   
        case 1 % Record Trace Bin
            [t2trace, t2time_vector, t2timestep, t2code] = readlabviewbinaries(strcat(handles.t2pathname,handles.t2filename));
            if get(handles.t2mult,'Value') == 1
                t2trace = handles.t2traceraw.*10; %%%% BAD
            end
            t2trace = t2trace.*1000;
        case 2 % ATF file
            fid = fopen(strcat(handles.t2pathname,handles.t2filename));
            if fid > 0 % opened
                for i = 1:10 % header to be skipped
                    line = fgetl(fid);
                end
                fffnoise = fscanf(fid, '%f\t%f', [2 inf]);
                
                t2timestep = fffnoise(1,3)-fffnoise(1,2);
                t2trace = fffnoise(2,:); % in pA
            end
            fclose(fid);
        case 3% Text file with two columns, same as above but without header
            fid = fopen(strcat(handles.t2pathname,handles.t2filename));
            if fid > 0
                fffnoise = fscanf(fid, '%f\t%f', [2 inf]);
                t2timestep = fffnoise(1,3)-fffnoise(1,2); %% column 1 is time
                t2trace = fffnoise(2,:).*1000; %% column 2 is current
            end
            fclose(fid);
        case 4 % ABF file
            segment_size_in_sec = str2num(get(handles.text_t1_duration,'String'));
            start_time = str2num(get(handles.text_t1_start,'String'));
                        
            [t2trace,si,lActualAcqLength]=abfload(strcat(handles.t2pathname,handles.t2filename),'start',start_time,'stop',segment_size_in_sec);
            
            t2timestep = si*10^(-6); % convert to sec
            t2trace = transpose(t2trace); % keep in pA, and switch to correct form
        otherwise
    end
    handles.t2timestep = t2timestep;
    handles.t2traceraw = t2trace;
    set(handles.t2time,'String',strcat(num2str(handles.t2timestep*length(handles.t2traceraw)),' seconds long'))
    set(handles.statustext,'String','Trace 2 File Loaded')
    %t1refresh(handles);
    guidata(gcbo, handles); % update handles


function t2fname_Callback(hObject, eventdata, handles)
% hObject    handle to t2fname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t2fname as text
%        str2double(get(hObject,'String')) returns contents of t2fname as a double


% --- Executes during object creation, after setting all properties.
function t2fname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t2fname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in t2plot.
function t2plot_Callback(hObject, eventdata, handles)
% hObject    handle to t2plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    t2trace = handles.t2traceraw;
    if get(handles.t2mult,'Value') == 1
        t2trace = handles.t2traceraw.*10;
    end
    if get(handles.t2dc,'Value') == 1
        t2trace = t2trace-mean(t2trace);
    end
    FsG = 1/(handles.t2timestep);
    [fGn PxxGn] = calcPSD(FsG, t2trace);
    if get(handles.t2smooth,'Value') == 1
        [smoothed_fGn smoothed_PxxGn] = smoothing(fGn, PxxGn);
        fGn = smoothed_fGn; 
        PxxGn = smoothed_PxxGn;
    end
    handles.t2fGn = fGn;
    handles.t2PxxGn = PxxGn;
    t2plotty(handles, 'here')
    guidata(hObject, handles);
    
function t2plotty(handles, location)
    if get(handles.lockaxes,'Value') == 1
        lxmin = str2num(get(handles.xmin,'String'));
        lxmax = str2num(get(handles.xmax,'String'));
        lymin = str2num(get(handles.ymin,'String'));
        lymax = str2num(get(handles.ymax,'String'));
    end
    if strcmp(location, 'here')
        cla(handles.axes1,'reset')
        axes(handles.axes1)
    end
    if strcmp(location, 'new')
        figure()
    end
    loglog(handles.t2fGn,handles.t2PxxGn,'-g')
    xlabel('Frequency (Hz)')
    ylabel('PSD (nA^2/Hz)')  
    if get(handles.lockaxes,'Value') == 1
        axis([lxmin lxmax lymin lymax]);
    end
    xlimits = get(handles.axes1,'XLim');
    ylimits = get(handles.axes1,'YLim');
    
    set(handles.xmin,'String',num2str(xlimits(1)))
    set(handles.xmax,'String',num2str(xlimits(2)))
    set(handles.ymin,'String',num2str(ylimits(1)))
    set(handles.ymax,'String',num2str(ylimits(2)))

% --- Executes on button press in t2enable.
function t2enable_Callback(hObject, eventdata, handles)
% hObject    handle to t2enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of t2enable


% --- Executes on selection change in t2pop.
function t2pop_Callback(hObject, eventdata, handles)
% hObject    handle to t2pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns t2pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from t2pop


% --- Executes during object creation, after setting all properties.
function t2pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t2pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in t2smooth.
function t2smooth_Callback(hObject, eventdata, handles)
% hObject    handle to t2smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of t2smooth


% --- Executes on button press in t3mult.
function t3mult_Callback(hObject, eventdata, handles)
% hObject    handle to t3mult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of t3mult


% --- Executes on button press in t3smooth.
function t3smooth_Callback(hObject, eventdata, handles)
% hObject    handle to t3smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of t3smooth


% --- Executes on button press in t3dc.
function t3dc_Callback(hObject, eventdata, handles)
% hObject    handle to t3dc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of t3dc


% --- Executes on button press in t3load.
function t3load_Callback(hObject, eventdata, handles)
% hObject    handle to t3load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [handles.t3filename, handles.t3pathname] = uigetfile( ...
    {  '*',  'All Files (*.*)'}, ...
       'Pick a file');
    set(handles.t3fname,'String',handles.t3filename)
    switch get(handles.t3pop,'Value')   
        case 1 % Record Trace Bin
            [t3trace, t3time_vector, t3timestep, t3code] = readlabviewbinaries(strcat(handles.t3pathname,handles.t3filename));
            if get(handles.t3mult,'Value') == 1
                t3trace = handles.t3traceraw.*10;
            end
        case 2 % ATF file
            fid = fopen(strcat(handles.t3pathname,handles.t3filename));
            if fid > 0 % opened
                for i = 1:10 % header to be skipped
                    line = fgetl(fid);
                end
                fffnoise = fscanf(fid, '%f\t%f', [2 inf]);
                
                t3timestep = fffnoise(1,3)-fffnoise(1,2);
                t3trace = fffnoise(2,:)./1000; % in nA
            end
            fclose(fid);
        case 3 % Text file with two columns, same as above but without header
            fid = fopen(strcat(handles.t3pathname,handles.t3filename));
            if fid > 0
                fffnoise = fscanf(fid, '%f\t%f', [2 inf]);
                t3timestep = fffnoise(1,3)-fffnoise(1,2); %% column 1 is time
                t3trace = fffnoise(2,:); %% column 2 is current
            end
            fclose(fid);
        case 4 % ABF file
            segment_size_in_sec = str2num(get(handles.text_t1_duration,'String'));
            start_time = str2num(get(handles.text_t1_start,'String'));
                        
            [t3trace,si,lActualAcqLength]=abfload(strcat(handles.t3pathname,handles.t3filename),'start',start_time,'stop',segment_size_in_sec);
            
            t3timestep = si*10^(-6); % convert to sec
            t3trace = transpose(t3trace/1000); % convert to nA, and switch to correct form
        otherwise
    end
    handles.t3timestep = t3timestep;
    handles.t3traceraw = t3trace;
    set(handles.t3time,'String',strcat(num2str(handles.t3timestep*length(handles.t3traceraw)),' seconds long'))
    set(handles.statustext,'String','Trace 3 File Loaded')
    %t1refresh(handles);
    guidata(hObject, handles);


function t3fname_Callback(hObject, eventdata, handles)
% hObject    handle to t3fname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t3fname as text
%        str2double(get(hObject,'String')) returns contents of t3fname as a double


% --- Executes during object creation, after setting all properties.
function t3fname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t3fname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in t3plot.
function t3plot_Callback(hObject, eventdata, handles)
% hObject    handle to t3plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    t3trace = handles.t3traceraw;
    if get(handles.t3mult,'Value') == 1
        t3trace = handles.t3traceraw.*10;
    end
    if get(handles.t3dc,'Value') == 1
        t3trace = t3trace-mean(t3trace);
    end
    FsG = 1/(handles.t3timestep);
    [fGn PxxGn] = calcPSD(FsG, t3trace);
    if get(handles.t3smooth,'Value') == 1
        [smoothed_fGn smoothed_PxxGn] = smoothing(fGn, PxxGn);
        fGn = smoothed_fGn; 
        PxxGn = smoothed_PxxGn;
    end
    handles.t3fGn = fGn;
    handles.t3PxxGn = PxxGn;
    t3plotty(handles, 'here')
    guidata(hObject, handles);
    
function t3plotty(handles, location)
    if get(handles.lockaxes,'Value') == 1
        lxmin = str2num(get(handles.xmin,'String'));
        lxmax = str2num(get(handles.xmax,'String'));
        lymin = str2num(get(handles.ymin,'String'));
        lymax = str2num(get(handles.ymax,'String'));
    end
    if strcmp(location, 'here')
        cla(handles.axes1,'reset')
        axes(handles.axes1)
    end
    if strcmp(location, 'new')
        figure()
    end
    loglog(handles.t3fGn,handles.t3PxxGn,'-b')
    xlabel('Frequency (Hz)')
    ylabel('PSD (nA^2/Hz)')  
    if get(handles.lockaxes,'Value') == 1
        axis([lxmin lxmax lymin lymax]);
    end
    xlimits = get(handles.axes1,'XLim');
    ylimits = get(handles.axes1,'YLim');
    
    set(handles.xmin,'String',num2str(xlimits(1)))
    set(handles.xmax,'String',num2str(xlimits(2)))
    set(handles.ymin,'String',num2str(ylimits(1)))
    set(handles.ymax,'String',num2str(ylimits(2)))

% --- Executes on button press in t3enable.
function t3enable_Callback(hObject, eventdata, handles)
% hObject    handle to t3enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of t3enable


% --- Executes on selection change in t3pop.
function t3pop_Callback(hObject, eventdata, handles)
% hObject    handle to t3pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns t3pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from t3pop


% --- Executes during object creation, after setting all properties.
function t3pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t3pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in buttonplotall.
function buttonplotall_Callback(hObject, eventdata, handles)
% hObject    handle to buttonplotall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if get(handles.lockaxes,'Value') == 1
        lxmin = str2num(get(handles.xmin,'String'));
        lxmax = str2num(get(handles.xmax,'String'));
        lymin = str2num(get(handles.ymin,'String'));
        lymax = str2num(get(handles.ymax,'String'));
    end
    cla(handles.axes1,'reset')
    axes(handles.axes1)    
    if get(handles.t1enable,'Value') == 1
        loglog(handles.t1fGn,handles.t1PxxGn,'-r')
        hold on;
    end
    if get(handles.t2enable,'Value') == 1
        loglog(handles.t2fGn,handles.t2PxxGn,'-g')
        hold on;
    end
    if get(handles.t3enable,'Value') == 1
        loglog(handles.t3fGn,handles.t3PxxGn,'-b')
    end
    
    
    xlabel('Frequency (Hz)')
    ylabel('PSD (nA^2/Hz)')  
    if get(handles.lockaxes,'Value') == 1
        axis([lxmin lxmax lymin lymax]);
    end
    xlimits = get(handles.axes1,'XLim');
    ylimits = get(handles.axes1,'YLim');
    
    set(handles.xmin,'String',num2str(xlimits(1)))
    set(handles.xmax,'String',num2str(xlimits(2)))
    set(handles.ymin,'String',num2str(ylimits(1)))
    set(handles.ymax,'String',num2str(ylimits(2)))

    t1in = get(handles.t1enable,'Value');
    t2in = get(handles.t2enable,'Value');
    t3in = get(handles.t3enable,'Value');
    tracecount = 0; 
    if ~t1in && ~t2in && ~t3in %% 1
        tracecount = 0;    
    elseif ~t1in && ~t2in && t3in %% 2
        stringy1 = get(handles.t3fname,'String');
        tracecount = 1; 
    elseif ~t1in && t2in && ~t3in %% 3
        stringy1 = get(handles.t2fname,'String');
        tracecount = 1; 
    elseif ~t1in && t2in && t3in %% 4
        stringy1 = handles.t2filename;
        stringy2 = get(handles.t3fname,'String');
        tracecount = 2; 
    elseif t1in && ~t2in && ~t3in %% 5
        stringy1 = get(handles.t1fname,'String');
        tracecount = 1; 
    elseif t1in && ~t2in && t3in %% 6
        stringy1 = get(handles.t1fname,'String');
        stringy2 = get(handles.t3fname,'String');
        tracecount = 2; 
    elseif t1in && t2in && ~t3in %% 7
        stringy1 = get(handles.t1fname,'String');
        stringy2 = get(handles.t2fname,'String');
        tracecount = 2; 
    elseif t1in && t2in && t3in %% 8
        stringy1 = get(handles.t1fname,'String');
        stringy2 = get(handles.t2fname,'String');
        stringy3 = get(handles.t3fname,'String');
        tracecount = 3; 
    end
    
    if tracecount == 1
        legend(stringy1,'Location','SouthWest')
    elseif tracecount == 2
        legend(stringy1,stringy2,'Location','SouthWest')
    elseif tracecount == 3
        legend(stringy1,stringy2,stringy3,'Location','SouthWest')
    end


function statustext_Callback(hObject, eventdata, handles)
% hObject    handle to statustext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of statustext as text
%        str2double(get(hObject,'String')) returns contents of statustext as a double


% --- Executes during object creation, after setting all properties.
function statustext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statustext (see GCBO)
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
    figure()   
    if get(handles.lockaxes,'Value') == 1
        lxmin = str2num(get(handles.xmin,'String'));
        lxmax = str2num(get(handles.xmax,'String'));
        lymin = str2num(get(handles.ymin,'String'));
        lymax = str2num(get(handles.ymax,'String'));
    end
    if get(handles.t1enable,'Value') == 1
        loglog(handles.t1fGn,handles.t1PxxGn,'-r')
        hold on;
    end
    if get(handles.t2enable,'Value') == 1
        loglog(handles.t2fGn,handles.t2PxxGn,'-g')
        hold on;
    end
    if get(handles.t3enable,'Value') == 1
        loglog(handles.t3fGn,handles.t3PxxGn,'-b')
    end
    
    
    xlabel('Frequency (Hz)')
    ylabel('PSD (nA^2/Hz)')
    if get(handles.lockaxes,'Value') == 1
        axis([lxmin lxmax lymin lymax]);
    end
    
    t1in = get(handles.t1enable,'Value');
    t2in = get(handles.t2enable,'Value');
    t3in = get(handles.t3enable,'Value');
    tracecount = 0; 
    if ~t1in && ~t2in && ~t3in %% 1
        tracecount = 0;    
    elseif ~t1in && ~t2in && t3in %% 2
        stringy1 = get(handles.t3fname,'String');
        tracecount = 1; 
    elseif ~t1in && t2in && ~t3in %% 3
        stringy1 = get(handles.t2fname,'String');
        tracecount = 1; 
    elseif ~t1in && t2in && t3in %% 4
        stringy1 = handles.t2filename;
        stringy2 = get(handles.t3fname,'String');
        tracecount = 2; 
    elseif t1in && ~t2in && ~t3in %% 5
        stringy1 = get(handles.t1fname,'String');
        tracecount = 1; 
    elseif t1in && ~t2in && t3in %% 6
        stringy1 = get(handles.t1fname,'String');
        stringy2 = get(handles.t3fname,'String');
        tracecount = 2; 
    elseif t1in && t2in && ~t3in %% 7
        stringy1 = get(handles.t1fname,'String');
        stringy2 = get(handles.t2fname,'String');
        tracecount = 2; 
    elseif t1in && t2in && t3in %% 8
        stringy1 = get(handles.t1fname,'String');
        stringy2 = get(handles.t2fname,'String');
        stringy3 = get(handles.t3fname,'String');
        tracecount = 3; 
    end
    
    if tracecount == 1
        legend(stringy1,'Location','SouthWest')
    elseif tracecount == 2
        legend(stringy1,stringy2,'Location','SouthWest')
    elseif tracecount == 3
        legend(stringy1,stringy2,stringy3,'Location','SouthWest')
    end
    
    



% --- Executes on button press in t1refresh.
function t1refresh_Callback(hObject, eventdata, handles)
    % hObject    handle to t1refresh (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    t1refresh(handles)
    guidata(hObject, handles);



% --- Executes on button press in lockaxes.
function lockaxes_Callback(hObject, eventdata, handles)
% hObject    handle to lockaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lockaxes




% --- Executes on button press in button_dump.
function button_dump_Callback(hObject, eventdata, handles)
% hObject    handle to button_dump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'t1fGn')
    assignin('base', 't1fGn', handles.t1fGn)
end
if isfield(handles,'t1PxxGn')
    assignin('base', 't1PxxGn', handles.t1PxxGn)
end
if isfield(handles,'t2fGn')
    assignin('base', 't2fGn', handles.t2fGn)
end
if isfield(handles,'t2PxxGn')
    assignin('base', 't2PxxGn', handles.t2PxxGn)
end
if isfield(handles,'t3fGn')
    assignin('base', 't3fGn', handles.t3fGn)
end
if isfield(handles,'t3PxxGn')
    assignin('base', 't3PxxGn', handles.t3PxxGn)
end



function text_t3_duration_Callback(hObject, eventdata, handles)
% hObject    handle to text_t3_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_t3_duration as text
%        str2double(get(hObject,'String')) returns contents of text_t3_duration as a double


% --- Executes during object creation, after setting all properties.
function text_t3_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_t3_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_t3_start_Callback(hObject, eventdata, handles)
% hObject    handle to text_t3_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_t3_start as text
%        str2double(get(hObject,'String')) returns contents of text_t3_start as a double


% --- Executes during object creation, after setting all properties.
function text_t3_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_t3_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_t2_duration_Callback(hObject, eventdata, handles)
% hObject    handle to text_t2_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_t2_duration as text
%        str2double(get(hObject,'String')) returns contents of text_t2_duration as a double


% --- Executes during object creation, after setting all properties.
function text_t2_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_t2_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_t2_start_Callback(hObject, eventdata, handles)
% hObject    handle to text_t2_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_t2_start as text
%        str2double(get(hObject,'String')) returns contents of text_t2_start as a double


% --- Executes during object creation, after setting all properties.
function text_t2_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_t2_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_t1_duration_Callback(hObject, eventdata, handles)
% hObject    handle to text_t1_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_t1_duration as text
%        str2double(get(hObject,'String')) returns contents of text_t1_duration as a double


% --- Executes during object creation, after setting all properties.
function text_t1_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_t1_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function text_t1_start_Callback(hObject, eventdata, handles)
% hObject    handle to text_t1_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_t1_start as text
%        str2double(get(hObject,'String')) returns contents of text_t1_start as a double


% --- Executes during object creation, after setting all properties.
function text_t1_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_t1_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

