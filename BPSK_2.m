% BER vs. SNR for BPSK
clear all;close all;clc;
Ns = 1e5;% no. transmitted symbols
EbN0dB_sim = 0:12; %SNR (dB)
EbN0dB_ana = 0:0.01:12; %SNR for analysis
BPSK_sig_set = [-1 1]; %BPSK signal set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EbN0_sim = 10.^(EbN0dB_sim/10); %SNR for simulation
N0_sim = mean(BPSK_sig_set.^2)./EbN0_sim; % N0
b = round(rand(1, Ns)); % info bits
s = BPSK_sig_set(b+1); %transmitted signal points
for j = 1:length(N0_sim)
   w = sqrt(N0_sim(j)/2)*randn(1, Ns); %AWGN
   r = s + w; % received signal values
   b2 = zeros(1, Ns); % received bits
   one_ind = find(r > 0); % decision  of bit 1
   b2(one_ind) = 1;
   BER_sim(j) = length(find(b ~= b2))/Ns; % simulate BER
end

plot(EbN0dB_sim, log10(BER_sim), 'ro'); % simulated BERs
ylim([-5 0]); xlabel('E_b/N_0 (dB)'); ylabel('log_{10}BER');
hold on; grid on;

EbN0_ana = 10.^(EbN0dB_ana/10) %SNR for analysis
Q = @(x)0.5*erfc(x/sqrt(2)); %definition of Q function
BER_ana = Q(sqrt(2*EbN0_ana));
plot(EbN0dB_ana, log10(BER_ana), 'b'); % analytical BERs
