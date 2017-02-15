% Cyclic Redundancy Check - example
% RJP Nov 2011

clear all
% crc generation
%g = [1 0 1 1];    % Generator polynomial 
g=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 0 1 1 0 1 1];


nmk = length(g)-1;          % n-k  
k = 16;                     % nr of bits to transmit 
n = nmk+k;                  % code length


gr = g(end:-1:1);

%d=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
d = randi([0 1],1,k);

dap = [d, zeros(1,nmk)];                    % shift data by n-k to the left
dap=dap(end:-1:1);

[q, parity]  = gfdeconv(dap,gr);             % parity
dap(1:length(parity))=parity;

    % crc receiver
 
s=[parity d(end:-1:1)];
[q1,parity1]=gfdeconv(dap,gr);


% c  = dr;
% c(1:length(parity)) = parity;               % codeword  =[parity, info]     
% 
% %%%%% BSC with Pe %%%%%%
% error = +(rand(1,n) < Pe);                  % error sequence
% 
% r = rem(c + error,2);                       % received sequence
% 
% %%% Error detection i.e. calculate syndrome
% 
% [q, syndrome1] = gfdeconv(r,gr);
% 
% % compare with error free codeword
% [quot, syndrome2] = gfdeconv(c,gr);
% 
% 
