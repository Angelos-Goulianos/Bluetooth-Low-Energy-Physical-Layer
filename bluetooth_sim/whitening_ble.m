

%transmitter side scrambler

d_in=[0 1 0 1 1 0 1 0 1 0 1 1 0];
current_state=[1 0 1 0 1 1 1]; %initial state depending on BLE channel to be used.we now assume channel: 0x17=23
for i=1:length(d_in)
d_out(i)=xor(d_in(i),current_state(7));
next_state=circshift(current_state,[0 1]);
next_state(5)=xor(current_state(4),current_state(7));
current_state=next_state;
end

 %receiver side scrambler 
 
 current_state=[1 0 1 0 1 1 1];
 d_rin=d_out; % input to the receiver
 for i=1:length(d_in)
     d_rout(i)=xor(d_rin(i),current_state(7));
     next_state=circshift(current_state,[0 1]);
     next_state(5)=xor(current_state(4),current_state(7));
     current_state=next_state;
 end