%% Constellation diagram for QPSK
clear all;close all;clc;
Ns = 1e4;% no. transmitted symbols
EbN0dB_sim = 5; %SNR (dB)
QPSK_sig_set = [1+i -1+i 1-i -1-i]; %QPSK signal set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EbN0_sim = 10^(EbN0dB_sim/10); %SNR for simulation
N0_sim = mean(abs(QPSK_sig_set.^2)/2)/EbN0_sim; % N0
b = round(rand(1, 2*Ns)); % info bits
bp1 = b(1:2:length(b));
bp2 = b(2:2:length(b));
m = 2*bp1+bp2+1; %indices for QPSK signal points
s = QPSK_sig_set(m); %transmitted signal points
w = sqrt(N0_sim/2)*randn(1, Ns)+ i*sqrt(N0_sim/2)*randn(1, Ns); %AWGN
r = s + w; % received signal values
b2 = []; % received bits
for n = 1:Ns
    if real(r(n)) > 0 & imag(r(n)) > 0;
        b2 = [b2 0 0];
    elseif real(r(n)) <= 0 & imag(r(n)) > 0;
        b2 = [b2 0 1];
    elseif real(r(n)) > 0 & imag(r(n)) <= 0;
        b2 = [b2 1 0];
    else
        b2 = [b2 1 1];
    end
end

BER = length(find(b ~= 2))/(2*Ns);
plot(r,'+'); hold on; grid on;
plot(QPSK_sig_set,'ro');
xlabel('Re\{r_n\}'); ylabel('Im\{r_n}');
axis([-2.5 2.5 -2.5 2.5]);
text(-0.4, 2.2, strcat('Eb/N0 = ', num2str(EbN0dB_sim),'dB'));


