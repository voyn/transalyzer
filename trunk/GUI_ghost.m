function varargout = GUI_ghost(varargin)
% GUI_GHOST M-file for GUI_ghost.fig
%      GUI_GHOST, by itself, creates a new GUI_GHOST or raises the existing
%      singleton*.
%
%      H = GUI_GHOST returns the handle to a new GUI_GHOST or the handle to
%      the existing singleton*.
%
%      GUI_GHOST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_GHOST.M with the given input arguments.
%
%      GUI_GHOST('Property','Value',...) creates a new GUI_GHOST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_ghost_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_ghost_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_ghost

% Last Modified by GUIDE v2.5 18-May-2011 18:49:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_ghost_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_ghost_OutputFcn, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes just before GUI_ghost is made visible.
function GUI_ghost_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to GUI_ghost (see VARARGIN)

    % Choose default command line output for GUI_ghost
    handles.output = hObject;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isappdata(0,'HandleMainGUI')
        HandleMainGUI=getappdata(0,'HandleMainGUI');
        if isappdata(HandleMainGUI,'eventtrace')

            handles.eventtrace = getappdata(HandleMainGUI,'eventtrace');
            handles.timetemp = getappdata(HandleMainGUI,'timetemp');
            handles.eventtraceunfilt = getappdata(HandleMainGUI,'eventtraceunfilt');
            handles.extrapoints = getappdata(HandleMainGUI,'extrapoints');
            handles.eventlengthinpoints = getappdata(HandleMainGUI,'eventlengthinpoints');
            handles.timestep = getappdata(HandleMainGUI,'timestep');
            handles.fsample = getappdata(HandleMainGUI,'fsample');
            set(handles.edit_startfreq,'String',num2str(handles.fsample/2-1));
        end
        plotini(handles)
    end


    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes GUI_ghost wait for user response (see UIRESUME)
    % uiwait(handles.figure1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotini(handles)

    cla(handles.axesghost,'reset')
    axes(handles.axesghost)
    
    plot(handles.timetemp,handles.eventtraceunfilt,'-g')
    hold on
    plot(handles.timetemp,handles.eventtrace,'-b')
    xlabel('Time (ms)')
    ylabel('Current (nA)')
    scaleset = str2double(get(handles.edit_scale,'String'));
    minsig = (1-scaleset/100)*min(handles.eventtraceunfilt);
    maxsig = (1+scaleset/100)*max(handles.eventtraceunfilt);
    axis([0 max(handles.timetemp) minsig maxsig])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = GUI_ghost_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_exit_Callback(hObject, eventdata, handles)
    % hObject    handle to button_exit (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    delete(gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_ghost_Callback(hObject, eventdata, handles)
    % hObject    handle to button_ghost (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    filtersvect = fliplr(str2num(get(handles.edit_stopfreq,'String')):str2num(get(handles.edit_step,'String')):str2num(get(handles.edit_startfreq,'String')));
    %mov=avifile('traces\LPfilter','compression','Cinepak');
    cla(handles.axesghost,'reset')
    axes(handles.axesghost)
    filt_order = 5;
    filt_Wn = 50000;
    filt_type = 'low';
    pausetime = str2double(get(handles.edit_pause,'String'));
    scaleset = str2double(get(handles.edit_scale,'String'));
    Fs = handles.fsample;
    Nyquistfreq = Fs/2;
    reduced_freq = filtersvect./Nyquistfreq;
    set(handles.edit_status,'String','Ghosting...')
    maxtime = max(handles.timetemp);
    minsig = (1-scaleset/100)*min(handles.eventtraceunfilt);
    maxsig = (1+scaleset/100)*max(handles.eventtraceunfilt);
    for k = 1:length(filtersvect)
        cla(handles.axesghost,'reset')
        [filt_b,filt_a] = butter(filt_order,filtersvect(k)/Nyquistfreq,filt_type);
        %filt_signal_events_noise = filter(filt_b,filt_a,handles.eventtraceunfilt);
        filt_signal_events_noise = filtfilt(filt_b,filt_a,handles.eventtraceunfilt);

        %filt_signal_events_noise(1:handles.extrapoints/2) = ones(handles.extrapoints/2,1)*mean(handles.eventtrace(1:handles.extrapoints));
        heightvect(k) = min(filt_signal_events_noise);
        heightvectup(k) = max(filt_signal_events_noise);
        %minindex=find(filt_signal_events_noise((proteinlocation(1)-10):(proteinlocation(length(proteinlocation))+10))==min(filt_signal_events_noise((proteinlocation(1)-10):(proteinlocation(length(proteinlocation))+10))));
        minindex=find(filt_signal_events_noise==heightvect(k));
        maxindex=find(filt_signal_events_noise==heightvectup(k));
        timeheightvect(k) = handles.timetemp(minindex);
        timeheightvectup(k) = handles.timetemp(maxindex);
        plot(timeheightvect, heightvect, '.r');
        hold on;
        plot(timeheightvectup, heightvectup, '.c');
        hold on;
        plot(handles.timetemp, handles.eventtrace, '-b');
        hold on;

        plot(handles.timetemp, filt_signal_events_noise, '-g');
        xlabel('Time (sec)')
        ylabel('Current (nA)')
        axis([0 maxtime minsig maxsig])
         %maxtime = (sampling_rate*(proteinlocation(1,1)+20));
         %mintime = (sampling_rate*(proteinlocation(1,1)-20));
         %axis([(sampling_rate*(proteinlocation(1,1)-20)) (sampling_rate*(proteinlocation(1,1)+20)) 0.98*minsig maxsig])
        text((maxtime)*0.06,(maxsig-0.98*minsig)*0.05+0.98*minsig,['Fc = ',num2str(filtersvect(k))],'HorizontalAlignment','left','Margin',10);
        %text((maxtime)*0.1,(maxsig-0.98*minsig)*0.05+0.98*minsig,['LP Butterworth',num2str(filt_order),'^{th} order'],'HorizontalAlignment','left','Margin',10);

        M(k) = getframe;
        %f2=getframe(gcf); % gets the gcf  
        %mov=addframe(mov,f2); % adds the frame into mov  

        %clf
        pause(pausetime)
    end
    %mov=close(mov); % closes the mov 
    fps = 12;
    repeatmovie = 1;
    %movie(M,repeatmovie,fps)
    
    if get(handles.check_savemov,'Value') == 1
        % make sure you use the right codec
        %movie2avi(M, 'myPeaks.avi', 'compression', 'Indeo5');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_startfreq_Callback(hObject, eventdata, handles)
    % hObject    handle to edit_startfreq (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of edit_startfreq as text
    %        str2double(get(hObject,'String')) returns contents of edit_startfreq as a double


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_startfreq_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to edit_startfreq (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_stopfreq_Callback(hObject, eventdata, handles)
    % hObject    handle to edit_stopfreq (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of edit_stopfreq as text
    %        str2double(get(hObject,'String')) returns contents of edit_stopfreq as a double


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_stopfreq_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to edit_stopfreq (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_step_Callback(hObject, eventdata, handles)
    % hObject    handle to edit_step (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of edit_step as text
    %        str2double(get(hObject,'String')) returns contents of edit_step as a double


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_step_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to edit_step (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_status_Callback(hObject, eventdata, handles)
% hObject

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_status_CreateFcn(hObject, eventdata, handles)

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pause_Callback(hObject, eventdata, handles)
    % hObject    handle to edit_pause (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of edit_pause as text
    %        str2double(get(hObject,'String')) returns contents of edit_pause as a double


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_pause_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to edit_pause (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_savemov_Callback(hObject, eventdata, handles)
    % hObject    handle to check_savemov (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: get(hObject,'Value') returns toggle state of check_savemov


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit7_Callback(hObject, eventdata, handles)
    % hObject    handle to edit7 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of edit7 as text
    %        str2double(get(hObject,'String')) returns contents of edit7 as a double


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit7_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to edit7 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_scale_Callback(hObject, eventdata, handles)
    % hObject    handle to edit_scale (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of edit_scale as text
    %        str2double(get(hObject,'String')) returns contents of edit_scale as a double


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_scale_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to edit_scale (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


