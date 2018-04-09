clear all;close all;clc;
N = 16;% no. OFDM subcarriers
NCP = 4; % CP length
NOFDM = 1e3;
EbN0dB_sim = 0:12; %SNR (dB)
EbN0dB_ana = 0:0.01:12;
QPSK_sig_set = [1+i -1+i 1-i -1-i]; %QPSK signal set
h = 0.4.^(0:4); % discrete-time CIR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EbN0_sim = 10.^(EbN0dB_sim/10); %SNR for simulation
N0_sim = (mean(abs(QPSK_sig_set).^2)/2)./EbN0_sim; % N0
b = round(rand(1, 2*NOFDM*N)); % info bits
bp1 = b(1:2:length(b));
bp2 = b(2:2:length(b));
m = 2*bp1+bp2+1; %indices for QPSK signal points
S = QPSK_sig_set(m); %transmitted signal points
s = [];
for j = 1:NOFDM
   tmp = sqrt(N)*ifft(S((j-1)*N + 1:j*N));
   s = [s tmp(N-NCP+1:N) tmp];
end
tmp2 = conv(h,s);
for k = 1: length(N0_sim)
    w = sqrt(N0_sim(k)/2)*randn(1, NOFDM*(N+NCP))+ i*sqrt(N0_sim(k)/2)*randn(1, NOFDM*(N+NCP)); %AWGN
    r = tmp2(1:length(s)) + w; % received signal values
    H = 1/sqrt(N)*fft(h, N);
    R = [];
   for j = 1:NOFDM
      tmp = r((j-1)*(N+NCP)+NCP+1:j*(N+NCP));
      R = [ R 1/sqrt(N)*fft(tmp)./(sqrt(N)*H)];
   %s = [s tmp(N-NCP+1:N) tmp];
   end
   b2 = []; % received bits
   for n = 1:N*NOFDM
     if real(R(n)) > 0 & imag(R(n)) > 0;
        b2 = [b2 0 0];
     elseif real(R(n)) <= 0 & imag(R(n)) > 0;
        b2 = [b2 0 1];
     elseif real(R(n)) > 0 & imag(R(n)) <= 0;
        b2 = [b2 1 0];
     else
        b2 = [b2 1 1];
     end
   end
  BER_sim(k) = length(find(b ~= b2))/(2*N*NOFDM);
end;

plot(EbN0dB_sim, log10(BER_sim), 'ro'); % simulation BERs
ylim([-5 0]); xlabel('E_B/N_0 (dB)'); ylabel('log_{10}BER');
hold on; grid on;

EbN0_ana = 10.^(EbN0dB_ana/10); % SNR for analysis
Q = @(x) 0.5*erfc(x/sqrt(2)); % definition of Q function
BER_ana = 0;
for j = 1:N;
    BER_ana = BER_ana + 1/N*Q(sqrt(2*N*abs(H(j))^2*EbN0_ana));
end
plot(EbN0dB_ana,log10(BER_ana),'b'); % analytical BERs
