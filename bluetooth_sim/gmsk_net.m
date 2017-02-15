clear all;
nrz_data = [1 -1 1 -1 1 -1 1 -1]; %sample data
%nrz_data = randsrc(1,10); % produces random -1's and 1's.
Tb = 1; % bit duration
BT = 0.3; % BT product of filter 
sps = 32; % samples per symbol
Ts = Tb/sps; % sample period
t=(-2*Tb:Ts:2*Tb);
alpha = 2*pi*BT/(sqrt(log(2)));
% Q-function defined in Q.m file
gauss = 0.5*erfc(alpha*(2*t-0.5)/sqrt(2)) - 0.5*erfc(alpha*(2*t+0.5)/sqrt(2)); % impulse response of Gaussian filter
K=pi/2/sum(gauss); % normalize filter. ensure phase transitions of pi/2
gauss = K*gauss;
nrz = upsample(nrz_data,sps);
nrz_gauss = conv(gauss,nrz); % filter the nrz data
subplot(2,1,1);stem(nrz_data);title('NRZ bits');xlabel('Time');ylabel('Amplitude');
subplot(2,1,2);plot(nrz_gauss);title('NRZ bits, after Gaussian Filtering');xlabel('Time');ylabel('Amplitude');
figure
nrz_gauss1 = cumsum(nrz_gauss); % integrate the data.
nrz_gauss2 = exp(j*nrz_gauss1);
plot(imag(nrz_gauss2));title('I and Q channels of modulated NRZ');xlabel('Time');ylabel('Amplitude');hold on
plot(real(nrz_gauss2),'r');
tt=0:Ts:(length(nrz_gauss1)-1)*Ts;
Ibranch=cos(real(nrz_gauss2)).*sin(2*pi*2.4e9*tt);
Qbranch=sin(imag(nrz_gauss2)).*cos(2*pi*2.4e9*tt);
r=Ibranch-Qbranch;
% first test without noise , noise = Nil
%nrz_gauss2 = downsample(nrz_gauss2,2*sps);
noisy_real = real(nrz_gauss2);
noisy_imag = imag(nrz_gauss2);
%pass the I and Q through a Matched Filter
filt_noisy_real = matched_filter(noisy_real,Tb,sps);
filt_noisy_imag = matched_filter(noisy_imag,Tb,sps);
% obtaining the phase of the analog signal
phase = atan2(filt_noisy_imag,filt_noisy_real);
derivative = diff(phase);
figure
stem(derivative);
%----------------------------
%definition of the Matched-Filter
function ans = matched_filter(x,T,Samples);
    t = (-T:T/Samples:T); 
    BT = 0.3;
    alpha = 2*pi*BT/(sqrt(log(2)));
    Mfil = Q(alpha*(t-0.5)) - Q(alpha*(t+0.5)); % impulse response of Matched filter
    
    % need to scale the filter, so that there is a phase change of pi/2 for
    % every bit change.
    
    K = pi/2/sum(Mfil);
    Mfil = K*Mfil;
    
    %convolve the filter with the signal(I or Q);
    
    ans = conv(Mfil,x);
    
end
%---------------------------------------------------------------
%definition of the Q-function
function ans = Q(x)
ans = 0.5*erfc(x/sqrt(2));
end