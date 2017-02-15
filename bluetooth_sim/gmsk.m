
% Pulse duration, sampling frequency and Gausssian filtering

T=10^-6; % Pulse duration
BT=0.5; % Normalized Bandwidth of Gaussian shaped filter
fs=8/T;% Sampling frequency (8 samples /bit)
N=fs*T; %Number of samples/pulse
Ts=1/fs;% sampling interval
t=-T:Ts:T;
a=sqrt(log(2)/2)*T/BT;
g=sqrt(pi)/a*exp(-(-pi/a*t).^2); % Gaussian filter

% Convolution of data with Gaussian filter
data=[-1 -1 1 1 1 -1]; % -1 -1 1 -1 1 -1 -1 1 1 -1 1 1 -1 -1 1 -1 1 -1 -1]; 
data_bits=kron(data,ones(1,N)); % sample data_bits
shaped_data=conv(data_bits,g);

% Normalization so that max phase difference to a single '1' passing
% through  the filter is pi*h (pi/2) for GMSK
h=0.5; % for GMSK
fm=h/(2*T); % frequency deviation
dev_fr=2*pi*fm; % cyclic frequency deviation (pi/2T)
data_one=ones(1,N);
shaped_one=conv(data_one,g);
c=N/sum(shaped_one);
shaped_data=shaped_data*c;

                                                  % GMSK modulation
                                    
phase=zeros(1,length(shaped_data)+1); % phase=0 at the beginning
for i=2:length(phase)
    phase(i)=phase(i-1)+pi*h*shaped_data(i-1)/N;
end
Ibranch=cos(phase); 
Qbranch=sin(phase);
tx_waveform=Ibranch+1i*Qbranch;
                                                % GMSK demodulation
                                                
rx_phase=atan2(imag(tx_waveform),real(tx_waveform));

                                                   % Unwrap phase              
for i=2:length(rx_phase)
    difference=rx_phase(i)-rx_phase(i-1);
    if difference > pi
        rx_phase(i:end)=rx_phase(i:end)-2*pi;
    elseif difference < -pi
        rx_phase(i:end)=rx_phase(i:end)+2*pi;
    end
end
%  sampling_instances=2*N+1:N:length(phase)-N;
%  sample_values=rx_phase(sampling_instances);
%  sample_values=[0 sample_values];
%  derivative=diff(sample_values);
%  rx_bits=derivative>0;


                                        % Slope detector
       
time=0:Ts:(length(rx_phase)-1)*Ts;
first=rx_phase(1:2*N+1);
time1=time(1:length(first));
p=polyfit(time1,first,1);
slope(1)=p(1);
sampl_inst=2*N+1;
for j=2:length(data)
    time2=time(sampl_inst:sampl_inst+8);
    f=rx_phase(sampl_inst:sampl_inst+8);
    sampl_inst=sampl_inst+8;
    pr=polyfit(time2,f,1);
    slope(j)=pr(1);
    clear pr time2
end
rx_bits=slope>0;
                                    % GMSK plot
                 
 
% sample_values=tx_waveform(sampling_time_instances)

% %  GMSK transmitted waveform
% 
% clear
% fc=2.4e9;
% T=10^-6; % Pulse duration
% BT=0.5; % Normalized Bandwidth of Gaussian shaped filter
% Ts=1/(16*fc);
% fs=1/Ts;
% N=fs*T; %Number of samples/pulse
% t=-T:Ts:T;
% a=sqrt(log(2)/2)*T/BT;
% g=sqrt(pi)/a*exp(-(-pi/a*t).^2);
% data_bits=[1 1 -1 1 1 -1 -1 1 -1 1 -1 -1 1 1 -1 1 1 -1 -1 1 -1 1 -1 -1]; 
% data_bits=kron(data_bits,ones(1,N)); % sample data_bits
% shaped_data=conv(data_bits,g);
% h=1/2; % for GMSK
% fm=h/(2*T); % frequency deviation
% dev_fr=2*pi*fm; % cyclic frequency deviation (pi/2T)
% data_one=ones(1,N);
% shaped_one=conv(data_one,g);
% c=N/sum(shaped_one);
% shaped_data=shaped_data*c;
% phase=zeros(1,length(shaped_data)+1); % phase=0 at the beginning
% for i=2:length(phase)
%     phase(i)=phase(i-1)+pi*h*shaped_data(i-1)/N;
% end
% Ibranch=cos(phase); 
% Qbranch=sin(phase);
% tt=0:Ts:Ts*(length(phase)-1);
% tx_waveform1=Ibranch.*cos(2*pi*fc*tt);
% tx_waveform2=Qbranch.*sin(2*pi*fc*tt);
% tx_waveform=tx_waveform1-tx_waveform2;




