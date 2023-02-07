% <<<<<<<<<<<<<<<<<<<< ASK Modulation and Demodulation >>>>>>>>>>>>>>>>>>>>

clc, clear all, close all;


% ******************* Digital/Binary input information ********************

x = randi(2, [1,100], 'int32') - 1; % Binary information as stream of bits (binary signal 0 or 1)
x1 = x;   
N = length(x);
Tb = 0.0001;   %Data rate = 10kHz i.e., bit period (second)
disp('Binary Input Information at Transmitter: ');
disp(x);


% ************* Represent input information as digital signal *************

nb = 20;   % Digital signal per bit
digit = []; 
for n = 1:1:N
    if x(n) == 1;
       sig = ones(1,nb);
    else x(n) == 0;
        sig = zeros(1,nb);
    end
     digit = [digit sig];
end

t1 = Tb/nb:Tb/nb:nb*N*(Tb/nb);   % Time period
figure('Name','ASK Modulation and Demodulation','NumberTitle','off');
subplot(4,1,1);
plot(t1,digit,'LineWidth',2.5);
grid on;
axis([0 Tb*N -0.5 1.5]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('Digital Input Signal');


% **************************** ASK Modulation *****************************

Ac1 = 1;     % Carrier amplitude for binary input '1'
Ac2 = -1;      % Carrier amplitude for binary input '0'
N0 = 5  ;
br = 1/Tb;    % Bit rate
Fc = br;   % Carrier frequency 
t2 = Tb/nb:Tb/nb:Tb;   % Signal time             

mod = [];
for (i = 1:1:N)
    if (x(i) == 1)
        y = Ac1*cos(2*pi*Fc*t2);   % Modulation signal with carrier signal 1
    else
        y = Ac2*cos(2*pi*Fc*t2);   % Modulation signal with carrier signal 2
    end
    mod = [mod y];
end

t3 = Tb/nb:Tb/nb:Tb*N;   % Time period
subplot(4,1,2);
plot(t3,mod);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('ASK Modulated Signal');


% ********************* Transmitted signal x ******************************

x = mod;


% ********************* Channel model h and w *****************************

h = 1;   % Signal fading 
w = sqrt(N0/2) * normrnd(0, 1, [1,N*nb]);   % Noise

% ********************* Received signal y *********************************

y = h.*x + w;   % Convolution
% y = awgn(x,1);
t3 = Tb/nb:Tb/nb:Tb*N;   % Time period
subplot(4,1,3);
plot(t3,y);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('ASK Modulated Signal');

% *************************** ASK Demodulation ****************************

s = length(t2);
demod = [];
for n = s:s:length(y)
  t4 = Tb/nb:Tb/nb:Tb;    % Time period
  c = cos(2*pi*Fc*t4)*sqrt(2/Tb);    % Carrier signal 
  mm = c.*y((n-(s-1)):n); % Convolution 
  t5 = Tb/nb:Tb/nb:Tb;
  z = trapz(t5,mm);       % Intregation 
  rz = round((sqrt(2/Tb)*z));
  Ac = ((Ac1 + Ac2)/2);   % Average of carrier amplitudes
  if(rz > Ac)             % Logical condition 
    a = 1;
  else
    a = 0;
  end
  demod = [demod a];
end

disp('Demodulated Binary Information at Receiver: ');
disp(demod);

% ********** Represent demodulated information as digital signal **********


digit = [];
for n = 1:length(demod);
    if demod(n) == 1;
       sig = ones(1,nb);
    else demod(n) == 0;
        sig = zeros(1,nb);
    end
     digit = [digit sig];
end

t5 = Tb/nb:Tb/nb:nb*length(demod)*(Tb/nb);   % Time period
subplot(3,1,3)
plot(t5,digit,'LineWidth',2.5);grid on;
axis([0 Tb*length(demod) -0.5 1.5]);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('ASK Demodulated Signal');
num = xor(x1, demod);
disp(sum(num));
% ************************** End of the program ***************************