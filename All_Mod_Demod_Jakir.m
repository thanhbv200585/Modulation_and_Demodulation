clc; clear all; close all;
carrier_frequency=5; %Hz
carrier_frequency2=7; %Hz This will be utilized for FSK's second carrier
x=[1 0 0 1 1 0 0 1] % input signal ;
length=size(x,2);
 
i=1;
while i<length+1
     t = i:0.001:i+1;
    if x(i)==1
       ask=sin(2*pi*carrier_frequency*t);
       fsk=sin(2*pi*carrier_frequency*t);
       psk=sin(2*pi*carrier_frequency*t);
    else
        ask=0;
        fsk=sin(2*pi*carrier_frequency2*t);
        psk=sin(2*pi*carrier_frequency*t+pi);
    end
   
    
    subplot(3,1,1);
    plot(t,ask);
    hold on;
    grid on;
    axis([1 10 -1 1]);
 title('Amplitude Shift Key')
    subplot(3,1,2);
    plot(t,fsk);
    hold on;
    grid on;
    axis([1 10 -1 1]);
title('Frequency Shift Key')
 
    subplot(3,1,3);
    plot(t,psk);
    hold on;
    grid on;
    axis([1 10 -1 1]);
 title('Phase Shift Key')
 
    i=i+1;
 end