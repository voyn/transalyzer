function varargout = GUI_gaussian(varargin)
% GUI_GAUSSIAN M-file for GUI_gaussian.fig
%      GUI_GAUSSIAN, by itself, creates a new GUI_GAUSSIAN or raises the existing
%      singleton*.
%
%      H = GUI_GAUSSIAN returns the handle to a new GUI_GAUSSIAN or the handle to
%      the existing singleton*.
%
%      GUI_GAUSSIAN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_GAUSSIAN.M with the given input arguments.
%
%      GUI_GAUSSIAN('Property','Value',...) creates a new GUI_GAUSSIAN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_gaussian_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_gaussian_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_gaussian

% Last Modified by GUIDE v2.5 01-Nov-2012 10:35:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_gaussian_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_gaussian_OutputFcn, ...
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


% --- Executes just before GUI_gaussian is made visible.
function GUI_gaussian_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_gaussian (see VARARGIN)

% Choose default command line output for GUI_gaussian
handles.output = hObject;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isappdata(0,'HandleMainGUI')
    HandleMainGUI=getappdata(0,'HandleMainGUI');
    if isappdata(HandleMainGUI,'dwellvector')
        handles.amplitude=getappdata(HandleMainGUI,'amplitude');
        handles.check_conductance=getappdata(HandleMainGUI,'check_conductance');
        handles.voltage=getappdata(HandleMainGUI,'voltage'); % in Volts
        handles.dwellvector=getappdata(HandleMainGUI,'dwellvector');
        handles.pathname=getappdata(HandleMainGUI,'pathname');
        if isappdata(HandleMainGUI,'newtrace')
            handles.newtrace=transpose(getappdata(HandleMainGUI,'newtrace'));
            
        end
        handles.dwellbins=getappdata(HandleMainGUI,'dwellbins');
        set(handles.edit_dwellbins,'String',num2str(handles.dwellbins))
        handles.cbbins=getappdata(HandleMainGUI,'cbbins');
        set(handles.edit_cbbins,'String',num2str(handles.cbbins))
        handles.currentbins=getappdata(HandleMainGUI,'currentbins');
        set(handles.edit_currentbins,'String',num2str(handles.currentbins))
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_gaussian wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_gaussian_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_dwellbins_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dwellbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dwellbins as text
%        str2double(get(hObject,'String')) returns contents of edit_dwellbins as a double


% --- Executes during object creation, after setting all properties.
function edit_dwellbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dwellbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_exit.
function button_exit_Callback(hObject, eventdata, handles)
% hObject    handle to button_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(gcf)

% --- Executes on selection change in pop_fit.
function pop_fit_Callback(hObject, eventdata, handles)
% hObject    handle to pop_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_fit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_fit


% --- Executes during object creation, after setting all properties.
function pop_fit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_numgaussdwell_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numgaussdwell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numgaussdwell as text
%        str2double(get(hObject,'String')) returns contents of edit_numgaussdwell as a double


% --- Executes during object creation, after setting all properties.
function edit_numgaussdwell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numgaussdwell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_openplot.
function button_openplot_Callback(~, ~, handles)
% hObject    handle to button_openplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fun_plot_it_all(handles, 'new')

% --- Executes on button press in button_replot.
function button_replot_Callback(~, ~, handles)
% hObject    handle to button_replot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fun_plot_it_all(handles, 'old')


function edit_mu1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mu1 as text
%        str2double(get(hObject,'String')) returns contents of edit_mu1 as a double


% --- Executes during object creation, after setting all properties.
function edit_mu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fun_plot_it_all(handles, location)
    selection = get(handles.uinput,'SelectedObject');
    switch get(selection,'Tag')
        case 'radio_block'
            y_values = handles.amplitude; 
            nbins = str2num(get(handles.edit_cbbins,'String'));
            numgauss = str2double(get(handles.edit_numgausscb,'String'));
        case 'radio_dwell'
            y_values = handles.dwellvector;
            nbins = str2num(get(handles.edit_dwellbins,'String'));
            numgauss = str2double(get(handles.edit_numgaussdwell,'String'));
        case 'radio_current'
            if handles.check_conductance == 1
                y_values = handles.newtrace./handles.voltage;
            elseif handles.check_conductance == 2
                y_values = handles.newtrace;
            else
                y_values = handles.newtrace;
            end
            
            nbins = str2num(get(handles.edit_currentbins,'String'));
            numgauss = str2double(get(handles.edit_numgausscurrent,'String'));
    end
    
    if strcmp(location, 'old')
        cla(handles.axes1,'reset')
        axes(handles.axes1)
    end
    if strcmp(location, 'new')
        %ylimits = get(handles.axeshist,'YLim');
        %xlimits = get(handles.axeshist,'XLim');
        figure()
    end
    
    
    
    [n,xbin]=hist(y_values,nbins);
    n_sum = sum(n);
    
    if strcmp(get(selection,'Tag'), 'radio_current')
        
        n = log(n);
    end
            
    switch get(handles.pop_fit,'Value') %% another copy of this code in the gaussian button func
        case 1 % histfit
            
            mr = nanmean(y_values); % Estimates the parameter, MU, of the normal distribution.
            sr = nanstd(y_values);  % Estimates the parameter, SIGMA, of the normal distribution.

            x=(-3*sr+mr:0.1*sr:3*sr+mr)';% Evenly spaced samples of the expected data range.

            hh = bar(xbin,n,1); % Plots the histogram.  No gap between bars.
            if strcmp(get(selection,'Tag'), 'radio_current')
                set(hh, 'FaceColor', 'b', 'EdgeColor', 'b');
            else
            	set(hh, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k');
            end
            ylimits = get(handles.axes1,'YLim');
            xlimits = get(handles.axes1,'XLim');
            hold on
            xd = get(hh,'Xdata'); % Gets the x-data of the bins.

            rangex = max(xd(:)) - min(xd(:)); % Finds the range of this data.
            binwidth = rangex/nbins;    % Finds the width of each bin.
            
            row = sum(~isnan(y_values));

            y = normpdf(x,mr,sr);  
            y = row*(y*binwidth);   % Normalization necessary to overplot the histogram.
            
            if strcmp(get(selection,'Tag'), 'radio_current')
                hh1 = plot(x,log(y),'r-','LineWidth',2);     % Plots density line over histogram.
            else
                hh1 = plot(x,y,'r-','LineWidth',2);     % Plots density line over histogram.
            end
            set(handles.edit_mu1,'String',num2str(mr))
            set(handles.edit_si1,'String',num2str(sr))
            set(handles.edit_ro1,'String',num2str(row))
            set(handles.edit_bi1,'String',num2str(binwidth))
            
            axis([xlimits(1) xlimits(2) 0 ylimits(2)]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 2 % gmdistribution.fit
            
            % see example at
            % http://jasonboulet.com/2011/11/10/gaussian-mixture-models/
            
            options = statset('Display','final');
            fit_result=gmdistribution.fit(y_values,numgauss,'Options',options); %fit mixed gaussian might take a few minutes on slow machines
            
            gaussPdf = pdf(fit_result,xbin');
            A_norm = sum(gaussPdf);
            gaussPdf = gaussPdf/A_norm;

            % separating N Gaussians

            for id_gauss = 1:numgauss
                mu_Value(id_gauss)    = fit_result.mu(id_gauss);
                sigma_Value(id_gauss) = sqrt(fit_result.Sigma(1,1,id_gauss));
                proportions(id_gauss) = n_sum*fit_result.PComponents(id_gauss);
                gaussPdfi(:,id_gauss) = proportions(id_gauss)*normpdf(xbin,mu_Value(id_gauss),sigma_Value(id_gauss))/A_norm;
            end

        
            % plot original data
            hh = bar(xbin,n,1);
            if strcmp(get(selection,'Tag'), 'radio_current')
                set(hh, 'FaceColor', 'b', 'EdgeColor', 'b');
            else
            	set(hh, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k');
            end
            
            hold on
            ylimits = get(handles.axes1,'YLim');
            xlimits = get(handles.axes1,'XLim');
            
            xd = get(hh,'Xdata'); % Gets the x-data of the bins.
            rangex = max(xd(:)) - min(xd(:));
            binwidth = rangex/nbins;    % Finds the width of each bin.
            
            % plot mixed GMM fit
            
            if strcmp(get(selection,'Tag'), 'radio_current')
                plot(xbin, log(n_sum.*gaussPdf), 'k', 'linewidth', 3);
            else
            	plot(xbin, n_sum.*gaussPdf, 'k', 'linewidth', 3);
            end

            colorstring = {'r' 'g' 'c' 'm' 'k' 'k' 'k' 'k' 'k' 'k' 'k' 'k' 'k' 'k' 'k' 'k' 'k' 'k'};

            % plot each distribution
            for distcounter = 1:numgauss
                
                if strcmp(get(selection,'Tag'), 'radio_current')
                    line(xbin, log(gaussPdfi(:,distcounter)), 'Color', colorstring{distcounter}, 'linewidth', 2);  
                else
                    line(xbin, gaussPdfi(:,distcounter), 'Color', colorstring{distcounter}, 'linewidth', 2);  
                end
            end

            set(handles.edit_mu1,'String',num2str(mu_Value(1)))
            set(handles.edit_si1,'String',num2str(sigma_Value(1)))
            set(handles.edit_ro1,'String',num2str(proportions(1)))
            set(handles.edit_bi1,'String',num2str(binwidth))
            if numgauss > 1
                set(handles.edit_mu2,'String',num2str(mu_Value(2)))
                set(handles.edit_si2,'String',num2str(sigma_Value(2)))
                set(handles.edit_ro2,'String',num2str(proportions(2)))
                set(handles.edit_bi2,'String',num2str(binwidth))
            end
            if numgauss > 2
                set(handles.edit_mu3,'String',num2str(mu_Value(3)))
                set(handles.edit_si3,'String',num2str(sigma_Value(3)))
                set(handles.edit_ro3,'String',num2str(proportions(3)))
                set(handles.edit_bi3,'String',num2str(binwidth))
            end
            axis([xlimits(1) xlimits(2) 0 ylimits(2)]);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 3 % manual
            mu1 = str2double(get(handles.edit_mu1,'String'));
            si1 = str2double(get(handles.edit_si1,'String'));
            ro1 = str2double(get(handles.edit_ro1,'String'));
            bi1 = str2double(get(handles.edit_bi1,'String'));
            if numgauss > 1
                mu2 = str2double(get(handles.edit_mu2,'String'));
                si2 = str2double(get(handles.edit_si2,'String'));
                ro2 = str2double(get(handles.edit_ro2,'String'));
                bi2 = str2double(get(handles.edit_bi2,'String'));
            end
            if numgauss > 2
                mu3 = str2double(get(handles.edit_mu3,'String'));
                si3 = str2double(get(handles.edit_si3,'String'));
                ro3 = str2double(get(handles.edit_ro3,'String'));
                bi3 = str2double(get(handles.edit_bi3,'String'));
            end
            
            
            mr = nanmean(y_values); % Estimates the parameter, MU, of the normal distribution.
            sr = nanstd(y_values);  % Estimates the parameter, SIGMA, of the normal distribution.
            x=(-3*sr+mr:0.1*sr:3*sr+mr)';% Evenly spaced samples of the expected data range.

            hh = bar(xbin,n,1); % Plots the histogram.  No gap between bars.
            if strcmp(get(selection,'Tag'), 'radio_current')
                set(hh, 'FaceColor', 'b', 'EdgeColor', 'b');
            else
            	set(hh, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k');
            end
            ylimits = get(handles.axes1,'YLim');
            xlimits = get(handles.axes1,'XLim');
            hold on
            y = normpdf(x,mu1,si1);  
            y = ro1*(y*bi1);   % Normalization necessary to overplot the histogram.
            if strcmp(get(selection,'Tag'), 'radio_current')
                hh1 = plot(x,log(y),'r-','LineWidth',2);     % Plots density line over histogram.
            else
                hh1 = plot(x,y,'r-','LineWidth',2);     % Plots density line over histogram.
            end
            
            if numgauss > 1
                y2 = normpdf(x,mu2,si2);
                y2 = ro2*(y2*bi2);   % Normalization necessary to overplot the histogram.
                %hh1 = plot(x,y2,'g-','LineWidth',2);     % Plots density line over histogram.
                if strcmp(get(selection,'Tag'), 'radio_current')
                    hh1 = plot(x,log(y2),'g-','LineWidth',2);     % Plots density line over histogram.
                else
                    hh1 = plot(x,y2,'g-','LineWidth',2);     % Plots density line over histogram.
                end
            end
            if numgauss > 2
                y3 = normpdf(x,mu3,si3);
                y3 = ro3*(y3*bi3);   % Normalization necessary to overplot the histogram.
                %hh1 = plot(x,y3,'c-','LineWidth',2);     % Plots density line over histogram.
                if strcmp(get(selection,'Tag'), 'radio_current')
                    hh1 = plot(x,log(y3),'c-','LineWidth',2);     % Plots density line over histogram.
                else
                    hh1 = plot(x,y3,'c-','LineWidth',2);     % Plots density line over histogram.
                end
            end
            axis([xlimits(1) xlimits(2) 0 ylimits(2)]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 4 % multi gauss GMM
            
            % based on
            % http://www.aquaphoenix.com/lecture/matlab10/page2.html
            
            [em_thr,em_thr_behavior,P,mu_Value,sigma_Value,pdf_x,xx,pdf_xx,cdf_xx] = em_1dim(y_values, numgauss);

            hh = bar(xbin,n,1);
            if strcmp(get(selection,'Tag'), 'radio_current')
                set(hh, 'FaceColor', 'b', 'EdgeColor', 'b');
            else
            	set(hh, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k');
            end
            
            ylimits = get(handles.axes1,'YLim');
            xlimits = get(handles.axes1,'XLim');
            hold on
            xd = get(hh,'Xdata'); % Gets the x-data of the bins.
            rangex = max(xd(:)) - min(xd(:));
            binwidth = rangex/nbins;    % Finds the width of each bin.
            row = sum(~isnan(y_values))/numgauss;
            colorstring = {'-r' '-g' '-c' '-m' '-k' '-k' '-k' '-k' '-k' '-k' '-k' '-k' '-k' '-k' '-k' '-k' '-k' '-k'};
            for distcounter = 1:numgauss
                y = row*(normpdf(xx,mu_Value(distcounter),sigma_Value(distcounter))*binwidth);
                if strcmp(get(selection,'Tag'), 'radio_current')
                    plot(xx,log(y),colorstring{distcounter},'LineWidth',2);
                else
                    plot(xx,y,colorstring{distcounter},'LineWidth',2);
                end
                
            end
            set(handles.edit_mu1,'String',num2str(mu_Value(1)))
            set(handles.edit_si1,'String',num2str(sigma_Value(1)))
            set(handles.edit_ro1,'String',num2str(row))
            set(handles.edit_bi1,'String',num2str(binwidth))
            if numgauss > 1
                set(handles.edit_mu2,'String',num2str(mu_Value(2)))
                set(handles.edit_si2,'String',num2str(sigma_Value(2)))
                set(handles.edit_ro2,'String',num2str(row))
                set(handles.edit_bi2,'String',num2str(binwidth))
            end
            if numgauss > 2
                set(handles.edit_mu3,'String',num2str(mu_Value(3)))
                set(handles.edit_si3,'String',num2str(sigma_Value(3)))
                set(handles.edit_ro3,'String',num2str(row))
                set(handles.edit_bi3,'String',num2str(binwidth))
            end
            axis([xlimits(1) xlimits(2) 0 ylimits(2)]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
   switch get(selection,'Tag')
        case 'radio_block'
            if handles.check_conductance == 1
                xlabel('Conductance (nS)')
            elseif handles.check_conductance == 2
                xlabel('Relative Blockade')
            else
                xlabel('Current Blockade (nA)')
            end
            ylabel('Counts')
        case 'radio_dwell'
            xlabel('Dwell Time (ms)')
            ylabel('Counts')
        case 'radio_current'
            if handles.check_conductance == 1
                xlabel('Conductance (nS)')
            elseif handles.check_conductance == 2
                xlabel('Relative Blockade')
            else
                xlabel('Current Blockade (nA)')
            end
            ylabel('Log Counts')
    end

    
    



function edit_si1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_si1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_si1 as text
%        str2double(get(hObject,'String')) returns contents of edit_si1 as a double


% --- Executes during object creation, after setting all properties.
function edit_si1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_si1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ro1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ro1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ro1 as text
%        str2double(get(hObject,'String')) returns contents of edit_ro1 as a double


% --- Executes during object creation, after setting all properties.
function edit_ro1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ro1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bi1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bi1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bi1 as text
%        str2double(get(hObject,'String')) returns contents of edit_bi1 as a double


% --- Executes during object creation, after setting all properties.
function edit_bi1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bi1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_mu2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mu2 as text
%        str2double(get(hObject,'String')) returns contents of edit_mu2 as a double


% --- Executes during object creation, after setting all properties.
function edit_mu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_si2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_si2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_si2 as text
%        str2double(get(hObject,'String')) returns contents of edit_si2 as a double


% --- Executes during object creation, after setting all properties.
function edit_si2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_si2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ro2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ro2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ro2 as text
%        str2double(get(hObject,'String')) returns contents of edit_ro2 as a double


% --- Executes during object creation, after setting all properties.
function edit_ro2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ro2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bi2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bi2 as text
%        str2double(get(hObject,'String')) returns contents of edit_bi2 as a double


% --- Executes during object creation, after setting all properties.
function edit_bi2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_mu3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mu3 as text
%        str2double(get(hObject,'String')) returns contents of edit_mu3 as a double


% --- Executes during object creation, after setting all properties.
function edit_mu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_si3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_si3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_si3 as text
%        str2double(get(hObject,'String')) returns contents of edit_si3 as a double


% --- Executes during object creation, after setting all properties.
function edit_si3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_si3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ro3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ro3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ro3 as text
%        str2double(get(hObject,'String')) returns contents of edit_ro3 as a double


% --- Executes during object creation, after setting all properties.
function edit_ro3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ro3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bi3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bi3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bi3 as text
%        str2double(get(hObject,'String')) returns contents of edit_bi3 as a double


% --- Executes during object creation, after setting all properties.
function edit_bi3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bi3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in button_dumpvariables.
function button_dumpvariables_Callback(hObject, eventdata, handles)
% hObject    handle to button_dumpvariables (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if isfield(handles,'amplitude')
        assignin('base', 'amplitude_data', handles.amplitude)
    end
    if isfield(handles,'dwellvector')
        assignin('base', 'dwell_data', handles.dwellvector)
    end
    if isfield(handles,'dwellbins')
        assignin('base', 'dwellbins', handles.dwellbins)
    end
    if isfield(handles,'cbbins')
        assignin('base', 'cbbins', handles.cbbins)
    end
    if isfield(handles,'currentbins')
        assignin('base', 'currentbins', handles.currentbins)
    end
    if isfield(handles,'newtrace')
        assignin('base', 'newtrace', handles.newtrace)
    end

    



function edit_cbbins_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cbbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cbbins as text
%        str2double(get(hObject,'String')) returns contents of edit_cbbins as a double


% --- Executes during object creation, after setting all properties.
function edit_cbbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cbbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_currentbins_Callback(hObject, eventdata, handles)
% hObject    handle to edit_currentbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_currentbins as text
%        str2double(get(hObject,'String')) returns contents of edit_currentbins as a double


% --- Executes during object creation, after setting all properties.
function edit_currentbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_currentbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_autosaveplot.
function button_autosaveplot_Callback(hObject, eventdata, handles)
% hObject    handle to button_autosaveplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    dir_str = strcat(handles.pathname,'plots', filesep); % in this sub directory
    if(~exist(dir_str, 'dir'))
        mkdir(dir_str);
    end
    
    % save plot
    newFig = figure;
    fun_plot_it_all(handles, 's')
    selection = get(handles.uinput,'SelectedObject');
    switch get(selection,'Tag')
        case 'radio_block'
            saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_amplitude_hist_fit.fig'), 'fig'); % save as .fig
            saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_amplitude_hist_fit_large.png'), 'png');
            set(newFig, 'PaperPositionMode', 'auto');
            saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_amplitude_hist_fit_small.png'), 'png');
            print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_amplitude_hist_fit.eps')) 
        case 'radio_dwell'
            saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_dwell_hist_fit.fig'), 'fig'); % save as .fig
            saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_dwell_hist_fit_large.png'), 'png');
            set(newFig, 'PaperPositionMode', 'auto');
            saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_dwell_hist_fit_small.png'), 'png');
            print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_dwell_hist_fit.eps')) 
        case 'radio_current'
            saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_current_hist_fit.fig'), 'fig'); % save as .fig
            saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_current_hist_fit_large.png'), 'png');
            set(newFig, 'PaperPositionMode', 'auto');
            saveas(newFig, strcat(handles.pathname,'plots', filesep,'plot_current_hist_fit_small.png'), 'png');
            print(newFig,'-depsc2',strcat(handles.pathname,'plots', filesep,'plot_current_hist_fit.eps')) 
    end
    
    close(newFig)





function edit_numgausscb_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numgausscb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numgausscb as text
%        str2double(get(hObject,'String')) returns contents of edit_numgausscb as a double


% --- Executes during object creation, after setting all properties.
function edit_numgausscb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numgausscb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_numgausscurrent_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numgausscurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numgausscurrent as text
%        str2double(get(hObject,'String')) returns contents of edit_numgausscurrent as a double


% --- Executes during object creation, after setting all properties.
function edit_numgausscurrent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numgausscurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
