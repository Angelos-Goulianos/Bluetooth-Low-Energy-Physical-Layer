
%This script implements an end to end  BLE transmitter-receiver under AWGN noise at 2.4 GHz
%considering the Gaussian pulse shaping as defined in the standard (i.e.
%BT=0.5) it also implements a differential phase detector and computes
%the Packet Error Probability (PER) and the Bit Error probability (BER) for
%various packet sizes, packet numbers, SNR values and modulation index
%h.According to standard h shoud be as follows: 0.45<h<0.55. The h=1/2
%value stands for GMSK modulation

%———————————————————————————


function [PER,BER,data,rx_packet]=ble_phy_pd(No_packets,No_bytes,SNR,h); % pd stands for :'phase-difference-based' receiver
% No_packets=1;
% No_bytes=1;
N_bits=No_bytes*8; % total bits 
data=rand(1,No_packets*N_bits)>0.5; % generate the data for all packets
packets=reshape(data,N_bits,No_packets); % each column is a packet
crc_pol=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 0 1 1 0 1 1]; % BLE CRC polynomial
crc_pol=crc_pol(end:-1:1); % write polynomial in reverse power due to MATLAB requirements for gfdecov function
zeros_length = length(crc_pol)-1; % length of zeros to add in the bits tx-ed
%SNR=10;
 for k=1:length(SNR)                                                                 
   f=SNR(k) % monitor progress of SNR loop
    
    for j=1:No_packets                                                                          
       packet=packets(:,j);
      
                                                              % CRC  insertion
                                                 
%      crc_form=[packet' zeros(1,zeros_length)];% crc packet starting from highest power
%      crc_form=crc_form(end:-1:1); % crc packet starting from lowest power 
%      [q, parity_tx] = gfdeconv(crc_form,crc_pol);
%      crc_packet=crc_form;
%      crc_packet(1:length(parity_tx))=parity_tx; % tx_ed packet with crc
%      [q1,parity1]=gfdeconv(crc_packet,crc_pol);
     
                                                              % Data   Whitening
                                                 
    current_state=[1 0 1 0 1 1 1]; %initial state depending on BLE channel to be used.we now assume channel: 0x17=23
    d_in=packet;
    for i=1:length(d_in)
        d_out(i)=xor(d_in(i),current_state(7));
        next_state=circshift(current_state,[0 1]);
        next_state(5)=xor(current_state(4),current_state(7));
        current_state=next_state;
    end     
    d_out=double(d_out);
    
                                                              % GFSK Modulation
                                                 
    T=10^-6; % Pulse duration (1/data_rate)
    BT=0.5; % Normalized Bandwidth of Gaussian shaped filter
    fs=8/T;% Sampling frequency 
    N=fs*T; %Number of samples/pulse
    Ts=1/fs;% sampling interval
    t=-T:Ts:T;
    a=sqrt(log(2)/2)*T/BT;
    g=sqrt(pi)/a*exp(-(-pi/a*t).^2); % Gaussian filter
    
    % Convolution of data with Gaussian filter
    d_out=2*d_out-1; % Form source coded NRZ data
    data_bits=kron(d_out,ones(1,N)); % sample data_bits
    shaped_data=conv(data_bits,g); % filter data with Gaussian pulse
    
    % Normalization so that [1/T x integral(-inf-inf)(data(one)*g(t))={(pi) x (h)}]
    %h=1/2;
    fm=h/(2*T); % frequency deviation
    dev_fr=2*pi*fm; % cyclic frequency deviation (pi/2T)
    data_one=ones(1,N);
    shaped_one=conv(data_one,g);
    c=N/sum(shaped_one);
    shaped_data=shaped_data*c;
    
    % Modulation
    phase=zeros(1,length(shaped_data)+1); % phase=0 at the beginning
    for i=2:length(phase)
        phase(i)=phase(i-1)+pi*h*shaped_data(i-1)/N;
    end
    Ibranch=cos(phase); 
    Qbranch=sin(phase);
    tx_waveform=Ibranch+1i*Qbranch;
    
                                                              % Additive white Gaussian Noise
    %ch = 1/sqrt(2)*[randn(1,length(tx_waveform)) + 1i*randn(1,length(tx_waveform))]; % white gaussian noise samples of 0dB variance 
    ch = 1/sqrt(2)*[randn(1,1) + 1i*randn(1,1)]; % white gaussian noise samples of 0dB variance 
    absol=abs(ch);
    tx_waveform = tx_waveform * ch;
                                                              
    rand('state', sum(100*clock));                                                                  
    S=10*log10(mean(abs(tx_waveform).^2)); % signal power  in dB                                         
    noise = 1/sqrt(2)*[randn(1,length(tx_waveform)) + 1i*randn(1,length(tx_waveform))]; % white gaussian noise samples of 0dB variance 
    noise_power=S-SNR(k); % noise variance-dB (noise power)
    noise_volt=10^(noise_power/20);
    tx_waveform=tx_waveform +noise_volt*noise; % noisy signal
                                         
                                                                  % GFSK demodulation
                                                
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
    
     %  Differential phase Detector
    sampling_instances=2*N+1:N:length(phase)-N;
    sample_values=rx_phase(sampling_instances);
    sample_values=[0 sample_values];% 0 is the reference phase
    derivative=diff(sample_values);
    rx_bits=derivative>0;% received ones and zeros
    
    
     % Slope detector
       
% time=0:Ts:(length(rx_phase)-1)*Ts;
% first=rx_phase(1:2*N+1);
% time1=time(1:length(first));
% p=polyfit(time1,first,1);
% slope(1)=p(1);
% sampl_inst=2*N+1;
% for j=2:length(d_out)
%     time2=time(sampl_inst:sampl_inst+8);
%     f=rx_phase(sampl_inst:sampl_inst+8);
%     sampl_inst=sampl_inst+8;
%     pr=polyfit(time2,f,1);
%     slope(j)=pr(1);
%     clear pr time2
% end
% rx_bits=slope>0;

%     derivative=diff(rx_phase);
%     rx_bits(1)=sum(derivative(1:2*N));
%     stop_point=2*N;
%     for jj=2:length(d_out)
%         rx_bits(jj)=sum(derivative(stop_point+1:stop_point+N));
%         stop_point=stop_point+N;
%     end
%     rx_bits=rx_bits>0;
 
%  derivative=diff(rx_phase);
%  sampling_instances=2*N:N:2*N+N*(length(d_out)-1);
%  sample_values=derivative(sampling_instances);
%  
%  rx_bits=sample_values>0;
        
  
                                                                   % Data De-Whitening
                                                 
    current_state=[1 0 1 0 1 1 1];
    d_rin=rx_bits; % input to the receiver
    for i=1:length(d_in)
        d_rout(i)=xor(d_rin(i),current_state(7));
        next_state=circshift(current_state,[0 1]);
        next_state(5)=xor(current_state(4),current_state(7));
        current_state=next_state;
    end
    
    drout=double(d_rout);% convert from logical to double
    
                                                                      % CRC receiver
                                                     
%     [q1,parity1]=gfdeconv(drout,crc_pol);
%     rx_data_packet=drout(end:-1:1);
 %   rx_packet=rx_data_packet(1:(length(packet)));
    rx_packet=drout;
    pe(j)=1-isequal(packet',rx_packet); %
    nErr(j)=size(find(packet'-rx_packet),2); % bits in error within a packet
%     equal=isequal(parity1,0);
%     equal=double(equal);
%     if equal==0 % if parity1 is not 0 then a packet in error
%         pe(j)=1;
%        elseif equal==1
%         pe(j)=0;
%        end  % end if
  % clear rx_bits derivative rx_phase rx_data_packet q1 drout parity1 rx_packet  equal 
    end % end packets
    
                                                                 % PER and BER calculation
                                                                 
    PER(k)=sum(pe)/No_packets; 
    BER(k)=sum(nErr)/(No_packets*N_bits);
                                                                              
 end % end SNR         
        
        
        