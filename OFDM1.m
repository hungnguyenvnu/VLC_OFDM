%% Signal values and constellation diagram for OFDM
clear all;close all;clc;
N = 16;% no. OFDM subcarriers
NCP = 4; % CP length
NOFDM = 1;
EbN0dB_sim = 100; %SNR (dB)
QPSK_sig_set = [1+i -1+i 1-i -1-i]; %QPSK signal set
h = 0.4.^(0:4); % discrete-time CIR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EbN0_sim = 10^(EbN0dB_sim/10); %SNR for simulation
N0_sim = (mean(abs(QPSK_sig_set).^2)/2)/EbN0_sim; % N0
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
tmp = conv(h,s);
w = sqrt(N0_sim/2)*randn(1, NOFDM*(N+NCP))+ i*sqrt(N0_sim/2)*randn(1, NOFDM*(N+NCP)); %AWGN
r = tmp(1:length(s)) + w; % received signal values
H = 1/sqrt(N)*fft(h, N);
R = [];
for j = 1:NOFDM
   tmp = r((j-1)*(N+NCP)+NCP+1:j*(N+NCP));
   R = [ R 1/sqrt(N)*fft(tmp)./(sqrt(N)*H)]
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

BER = length(find(b ~= b2))/(2*N*NOFDM);
k = 0:NOFDM*N-1;
k2 = 0:NOFDM*(N+NCP)-1;

figure(1); % transmit signal
subplot(2,2, 1), stem(k, real(S),'r');
xlabel('k'); ylabel('Re\{S_k\}'); xlim([0 NOFDM*N-1]);
subplot(2,2,2), stem(k, imag(S),'r');
xlabel('k'); ylabel('Im\{S_k\}'); xlim([0 NOFDM*N-1]);
subplot(2,2,3), stem(k2, real(s),'b');
xlabel('n'); ylabel('Re\{s_n\}'); xlim([0 NOFDM*(N+NCP)-1]);
subplot(2,2,4), stem(k2, imag(s),'b');
xlabel('n'); ylabel('Im\{s_n\}'); xlim([0 NOFDM*(N+NCP)-1]);

figure(2); % receive signal
subplot(2,2, 1), stem(k2, real(r),'b');
xlabel('n'); ylabel('Re\{r_k\}'); xlim([0 NOFDM*(N+NCP)-1]);
subplot(2,2,2), stem(k2, imag(r),'b');
xlabel('n'); ylabel('Im\{r_k\}'); xlim([0 NOFDM*(N+NCP)-1]);
subplot(2,2,3), stem(k, real(R),'r');
xlabel('k'); ylabel('Re\{R_k\}'); xlim([0 NOFDM*N-1]);
subplot(2,2,4), stem(k, imag(R),'r');
xlabel('k'); ylabel('Re\{R_k\}'); xlim([0 NOFDM*N-1]);

figure(3);
plot(R, '+'); hold on; grid on;
plot(QPSK_sig_set,'ro');
xlabel('Re\{R_k}'); ylabel('Im\{R_k}');
axis([-2.5 2.5 -2.5 2.5]);
text(-0.4, 2.2, strcat('Eb/N0=', num2str(EbN0dB_sim),'dB'))

