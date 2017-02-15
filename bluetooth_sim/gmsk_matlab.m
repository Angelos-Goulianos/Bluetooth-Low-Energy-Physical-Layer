
data=[0 0 1 1 1 zeros(1,16)]';
hMod = comm.GMSKModulator('BitInput',true, 'BandwidthTimeProduct',0.5);
modSignal = step(hMod, data);
hDemod = comm.GMSKDemodulator('BitOutput',true,'BandwidthTimeProduct',0.5);
demodSignal=step(hDemod,modSignal);
rx_signal=demodSignal(17:end);