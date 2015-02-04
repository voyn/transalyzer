function [f,Pxx] = calcPSD(Fs, XXX)

%first take the length of time vector
N = length(XXX);

%define the frequency of the power in L/2 units from 0 to freq/2 (niquist
%rate)
f =(Fs/2)*linspace(0,1,ceil(N/2));

%here calculate the power by taking the power of 2 of the absolute fft
% and divide it by L (amount of datapoints * frequency)

Pxx_temp = (fft(XXX).*conj(fft(XXX))/(N*Fs));

% throw away half the spectrum for the niquist rate

Pxx = 2*Pxx_temp(1:ceil(N/2));

Pars_time = norm(XXX).^2;
Pars_freq = norm(Pxx_temp).^2./numel(Pxx_temp)*Fs;
    
Ratio_Time_Freq = Pars_time/Pars_freq;

display(strcat('Parseval Ration Time / Freq = ',num2str(Ratio_Time_Freq)))

% %first i take the length in L
% L=length(D);
% Frequency=5*10^5;
% %define the frequency of the power in L/2 units from 0 to freq/2 (niquist
% %rate)
% freq2=(Frequency/2)*linspace(0,1,ceil(L/2));
% %here i calculate the power by taking the power of 2 of the absolute fft
% %and divide it by L (amount of datapoints * frequency)
% Ys=(  fft(D).*conj(fft(D))  /(length(D)*Frequency));
% %Ys I throw away half the spectrum for the niquist rate
% %I've thrown away half the spectrum (EXCEPT FOR THE 0 FREQUENCY), so I should multiply by 2 to get he
% %power
% %So never check your parseval relation with the power spectrum
% %The offset doesn't matter for the power spectrum since you do a logplot and offset ahs frequency zero
% %so that's why it doesn't matter
% Pl=2*Ys(1:ceil(L/2));
% %then i plot
% loglog(freq2,Pl);title('Power spectrum');xlabel('log Frequency');ylabel('log Power pA^2/hz');




%%%%%%%%%%%%%%% OLD BAD
% N = length(XXX);
% dt = 1/Fs;
% T = N*dt;
% df = 1/T;
% 
% f_temp     = df*(0:N-1)';   % max(f) = Fs - df 
% 
% Pxx_temp   = abs(fft(XXX)).^2/Fs/N;   % Power Spectral Density
% 
% Pxx = 2.*Pxx_temp(1:end/2+1);
% f = f_temp(1:end/2+1);
%%%%%%%%%%%%%%%

% nfft = 2^nextpow2(length(XXX));
% Pxx = abs(fft(XXX,nfft)).^2/length(XXX)/Fs;
% f = Fs*(0:nfft-1)/nfft;

% Create a single-sided spectrum
% Hpsd = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);  
% plot(Hpsd); 