function [trace, time_vector, timestep, code] = readchimlogfile(filename)

%%%%%%%%%%%%%%%%%%%% 
%
% Modified from CHIMERA_loganalysis
% function:  pushbutton_loadpreview_Callback(hObject, eventdata, handles)
% resampling and filtering cut out
% should be able to read in MATLAB .mat files to obtain settings
%
% DV 21/04/15
%
%%%%%%%%%%%%%%%%%%%%


% ~~~~~~~~~~~~
matfilename = strrep(filename,'.log','.mat');

% defaults
%SETUP_TIAgain=10e6;
%SETUP_preADCgain=0.97;
%SETUP_pAoffset=0;
%SETUP_mVoffset=0;
%SETUP_ADCVREF=2.5;
%SETUP_ADCBITS=14;
%ADCSAMPLERATE=4.1667e6;
% /defaults

if(exist(matfilename)==2)
   load(matfilename);
else
    msgbox('.mat file does not exist','Error','error')
end

samplerate = ADCSAMPLERATE;
TIAgain = SETUP_TIAgain;
preADCgain = SETUP_preADCgain;
currentoffset = SETUP_pAoffset;
voltageoffset = SETUP_mVoffset;
ADCvref = SETUP_ADCVREF;
ADCbits = SETUP_ADCBITS;

closedloop_gain = TIAgain*preADCgain;

% ~~~~~~~~~~~~

%LPfiltercutoff = 1e3*str2double(get(handles.edit_kHzbandwidth,'String'));
%%%%% Why resample data? %%%%%
%outputsamplerate = 1e3*str2double(get(handles.edit_outputsamplerate,'String'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ~~~~~~~~~~~~
%set(handles.text_status,'String','[LOADING]');
%set(handles.text_status,'BackgroundColor','red');
%plot(handles.axes_tracepreview,[0 1],[0 0],'r');
%xlabel('Time (s)');
%ylabel('Current (A)');
%pause(0.1);
% ~~~~~~~~~~~~


readfid = fopen(filename,'r');

%right-justified
%logdata = ( -ADCvref + (2*ADCvref) * double(mod(fread(readfid,'uint16'),2^ADCbits)) / 2^ADCbits ) / (TIAgain*preADCgain);

%left-justified
%     bitmask = (2^16 - 1) - (2^(16-ADCbits) - 1);
%     rawvalues = fread(readfid,'uint16');
%     readvalues = bitand(cast(rawvalues,'uint16'),bitmask);
%     logdata = ( -ADCvref + (2*ADCvref) * double(readvalues) / 2^16 ) / (TIAgain*preADCgain);
    
bitmask = (2^16 - 1) - (2^(16-ADCbits) - 1);
rawvalues = fread(readfid,'uint16');
readvalues = bitand(cast(rawvalues,'uint16'),bitmask);
%logdata = -ADCvref + (2*ADCvref) * double(readvalues) / 2^16;
logdata = ADCvref - (2*ADCvref) * double(readvalues) / 2^16;  % EDITED 11/25/13 to Match Current GUI
    
    
trace = (-logdata./closedloop_gain + currentoffset)'*1E9; % modified for Transalyzer

fclose(readfid);

outputsamplerate = samplerate;
timestep = 1/outputsamplerate;
time_vector = (1:length(trace))' ./ outputsamplerate;

% quick check and error message
if(length(trace)>1)
    code =1;
else
    code =0;
end

%filterorder = floor(samplerate/LPfiltercutoff*16);      % EDITED 8/15/2012
%myLPfilter = fir1(filterorder, LPfiltercutoff/(0.5*samplerate), 'low');
%if (get(handles.LPF_enable,'Value'))
%logdata = filter(myLPfilter,1,logdata);
%end

%%%% RESAMPLING DATA ???? %%%%%%%
%[P,Q] = rat(outputsamplerate/samplerate,0.02);
%set(handles.edit_outputsamplerate,'String',samplerate*P/Q*1e-3)
%%logdata = resample(logdata,P,Q);
%
%logdata = resample(logdata,P,Q,0);
%
%logdata = logdata(filterorder:(length(logdata)-filterorder));
%
%t = (1:length(logdata))' ./ outputsamplerate;
%
%plot(handles.axes_tracepreview,t*1e3,logdata*1e12);
%xlabel('Time (ms)');
%ylabel('Current (pA)');


end