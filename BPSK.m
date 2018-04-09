% Constellation diagram for BPSK
clear all;close all;clc;
Ns = 1e5;% no. transmitted symbols
EbN0dB_sim = 10; %SNR (dB)
BPSK_sig_set = [-1 1]; %BPSK signal set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EbN0_sim = 10^(EbN0dB_sim/10); %SNR for simulation
N0_sim = mean(BPSK_sig_set.^2)/EbN0_sim; % N0
b = round(rand(1, Ns)); % info bits
s = BPSK_sig_set(b+1); %transmitted signal points
w = sqrt(N0_sim/2)*randn(1, Ns); %AWGN
r = s + w; % received signal values
b2 = zeros(1, Ns); % received bits
one_ind = find(r > 0) % decision  of bit 1
b2(one_ind) = 1;
BER_sim = length(find(b ~= b2))/Ns % simulate BER

plot(r, zeros(1, Ns),'+'); hold on ; grid on;
plot(BPSK_sig_set, [0 0], 'ro');
xlabel('Re\{r_\n}')
ylabel('Im\{r_n}\')
axis([-2.5 2.5 -2.5 2.5]);
text(-0.4, 2.2, strcat('Eb/N0=', num2str(EbN0dB_sim),'dB'));