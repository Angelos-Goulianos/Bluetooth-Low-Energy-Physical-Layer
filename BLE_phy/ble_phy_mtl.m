
%This script implements an end to end  BLE transmitter-receiver under AWGN noise at 2.4 GHz
%considering the Gaussian pulse shaping as defined in the standard (i.e.
%BT=0.5) it also implements a differential phase detector and computes
%the Packet Error Probability (PER) and the Bit Error probability (BER) for
%various packet sizes, packet numbers, SNR values and modulation index
%h.According to standard h shoud be as follows: 0.45<h<0.55. The h=1/2
%value stands for GMSK modulation

function [PER,BER]=ble_phy_mtl(No_packets,No_bytes,SNR); % pd stands for :'phase-difference-based' receiver
% No_packets=1;
% No_bytes=1;
hMod = comm.GMSKModulator('BitInput',true, 'BandwidthTimeProduct',0.5);
hDemod = comm.GMSKDemodulator('BitOutput',true,'BandwidthTimeProduct',0.5);

N_bits=No_bytes*8; % bits in a packet
data=rand(1,No_packets*N_bits)>0.5; % generate the data for all packets
packets=reshape(data,N_bits,No_packets); % each column is a packet
crc_pol=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 0 1 1 0 1 1]; % BLE CRC polynomial
crc_pol=crc_pol(end:-1:1); % write polynomial in reverse power due to MATLAB requirements for gfdecov function
zeros_length = length(crc_pol)-1; % length of zeros to add in the bits tx-ed
%SNR=10;
 for k=1:length(SNR)                                                                 
   f=SNR(k) % mmonitor progress of SNR loop
    
    for j=1:No_packets                                                                          
       packet=packets(:,j);
      
                                                             % CRC  insertion
                                                 
     crc_form=[packet' zeros(1,zeros_length)];% crc packet starting from highest power
     crc_form=crc_form(end:-1:1); % crc packet starting from lowest power 
     [q, parity_tx] = gfdeconv(crc_form,crc_pol);
     crc_packet=crc_form;
     crc_packet(1:length(parity_tx))=parity_tx; % tx_ed packet with crc 
     %[q1,parity1]=gfdeconv(crc_packet,crc_pol);
     
                                                             % Data Whitening
                                                 
    current_state=[1 0 1 0 1 1 1]; %initial state depending on BLE channel to be used.we now assume channel: 0x17=23
    d_in=crc_packet;
    for i=1:length(d_in)
        d_out(i)=xor(d_in(i),current_state(7));
        next_state=circshift(current_state,[0 1]);
        next_state(5)=xor(current_state(4),current_state(7));
        current_state=next_state;
    end     
    d_out=double(d_out);
    
                                                            % GFSK Modulation
                                                                  
    data=[d_out zeros(1,16)]'; % will demodulate 16 bits less cause of the tracebacklength which is 16 and  the receiver starts decoding after the 17-th bit
    tx_waveform = step(hMod, data);
    
                                                       % Additive white Gaussian Noise
                                                                         
    rand('state', sum(100*clock));                                                                  
    S=10*log10(mean(abs(tx_waveform).^2)); % signal power  in dB                                         
    noise = 1/sqrt(2)*[randn(1,length(tx_waveform)) + 1i*randn(1,length(tx_waveform))]; % white gaussian noise samples of 0dB variance 
    noise_power=S-SNR(k); % noise variance-dB (noise power)
    noise_volt=10^(noise_power/20);
    tx_waveform=tx_waveform +noise_volt*noise.'; % noisy signal
                                         
                                                             % GFSK demodulation
                                                                  
    demodSignal=step(hDemod,tx_waveform);
    rx_bits=demodSignal(17:end);
                                                
    
                                                              % Data De-Whitening
                                                 
    current_state=[1 0 1 0 1 1 1];
    d_rin=rx_bits'; % input to the receiver
    for i=1:length(d_in)
        d_rout(i)=xor(d_rin(i),current_state(7));
        next_state=circshift(current_state,[0 1]);
        next_state(5)=xor(current_state(4),current_state(7));
        current_state=next_state;
    end
    
    drout=double(d_rout);% convert from logical to double
    
                                                               % CRC receiver
                                                     
    [q1,parity1]=gfdeconv(drout,crc_pol);
    rx_data_packet=drout(end:-1:1);
    rx_packet=rx_data_packet(1:(length(packet)));
    %pe(j)=1-isequal(packet',rx_packet); %
    nErr(j)=size(find(packet'-rx_packet),2); % bits in error within a packet
    equal=isequal(parity1,0);
    if equal==0 % if parity1 is not 0 then a packet in error
        pe(j)=1;
       elseif equal==1
        pe(j)=0;
       end  % end if
%    clear rx_bits derivative rx_phase rx_packet   
    end % end packets
    
                                                          % PER and BER calculation
                                                                 
    PER(k)=sum(pe)/No_packets; 
    BER(k)=sum(nErr)/(No_packets*N_bits);
                                                                              
 end % end SNR         
        
        
        