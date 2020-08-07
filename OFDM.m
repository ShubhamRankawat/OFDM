close all;
clear all;
rng('shuffle');

SNRdB = [1:5:45];
Nsub= 512;
Ncp = round (Nsub/10);
numBlocks = 10000;
numTaps = 2;
BER = zeros(size(SNRdB));
SNR = zeros(size(SNRdB));
for L = 1:numBlocks
bits = randi([0,1],[1,Nsub]);
 h=1/sqrt(2)*(randn(1,numTaps)+j*randn(1,numTaps));
   Hfreq = fft(h,Nsub);
   ChNoise = (randn(1,numTaps+ Nsub + Ncp -1)+j*randn(1,numTaps+ Nsub + Ncp -1));
for K=1:length(SNRdB)
    SNR(K)=10^(SNRdB(K)/10);
    Loadedbits = sqrt(SNR(K))*(2*bits-1);
    TxSamples = ifft(Loadedbits);
    TxSamplesCp = [TxSamples(Nsub-Ncp+1:Nsub),TxSamples];
    Rxbits = conv(h,TxSamplesCp)+ ChNoise;
    RxbitsWithoutCp = Rxbits(Ncp + 1:Ncp + Nsub);
    RxbitsFFT = fft(RxbitsWithoutCp,Nsub);
    ProcessedBits=RxbitsFFT./Hfreq;
    DecodedBits = ((real(ProcessedBits))>=0);
    BER(K)= BER(K) + sum(DecodedBits~=bits);
end
end
eSNR = numTaps*SNR/Nsub;
BER = BER/(numBlocks*Nsub);
semilogy(SNRdB,BER,'b s','linewidth',2.0);
hold on;
semilogy(SNRdB,0.5*(1-sqrt(eSNR./(2+eSNR))),'r-.','linewidth',2.0);
axis tight;
grid on;
legend('OFDM','Theory')
xlabel('SNR(dB)');
ylabel('BER');
title('BER vs SNR(dB)');

