%% BER and SNR for QPSK
clear all;close all;clc;
Ns = 1e4;% no. transmitted symbols
EbN0dB_sim = 0:12; %SNR (dB)
EbN0dB_ana = 0:0.01:12; % SNR for analysis
QPSK_sig_set = [1+i -1+i 1-i -1-i]; %QPSK signal set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EbN0_sim = 10.^(EbN0dB_sim/10); %SNR for simulation
N0_sim = (mean(abs(QPSK_sig_set).^2)/2)./EbN0_sim; % N0
b = round(rand(1, 2*Ns)); % info bits
bp1 = b(1:2:length(b));
bp2 = b(2:2:length(b));
m = 2*bp1+bp2+1; %indices for QPSK signal points
s = QPSK_sig_set(m); %transmitted signal points
for j = 1:length(N0_sim)
    w = sqrt(N0_sim(j)/2)*randn(1, Ns)+ i*sqrt(N0_sim(j)/2)*randn(1, Ns); %AWGN
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
    BER_sim(j) = length(find(b ~= b2))/(2*Ns);
end
plot(EbN0dB_sim, log10(BER_sim), 'ro'); % simulated BERs
ylim([-5 0]); xlabel('E_b/N_0 (dB)'); ylabel('log_{10}BER');
hold on; grid on;

EbN0_ana = 10.^(EbN0dB_ana/10) %SNR for analysis
Q = @(x)0.5*erfc(x/sqrt(2)); %definition of Q function
BER_ana = Q(sqrt(2*EbN0_ana));
plot(EbN0dB_ana, log10(BER_ana), 'b'); % analytical BERs



