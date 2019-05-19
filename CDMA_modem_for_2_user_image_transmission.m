function cdmamodem(user1,user2,snr_in_dbs)
% >>>multiple access b/w 2 users using DS CDMA
%% code by MAHIMA S. TENDULKAR
%snr_in_dbs=-50

close all
clearvars;
snr_in_dbs=50;

%code for image input1
inA= imread('prologo1low.jpg');
figure();
imshow(inA);
in11=imresize(inA,0.05);
N=numel(in11);
in1=reshape(in11,N,1);
bin1=de2bi(in1,'left-msb');
input1=reshape(bin1',numel(bin1),1);
len1=length(input1);
figure()
imshow(in11);
%code for image input2
inB= imread('prologo2low.jpg');
figure();
imshow(inB);
in22=imresize(inB,0.05);
M=numel(in22);
in2=reshape(in22,M,1);
bin2=de2bi(in2,'left-msb');
input2=reshape(bin2',numel(bin2),1);
len2=length(input2);
figure();
imshow(in22);

Z=max(len1,len2);

%Zero Padding

if len1<Z
    
    input1(Z,1)=0;
    
else if len2<Z
        
        input2(Z,1)=0;
    end 
end
inputA=double(input1);
inputB=double(input2);
user1=inputA;
user2=inputB;

%% to convert the binary sequences to bipolar NRZ format

length_user1=length(user1);
length_user2=length(user2);
for i=1:length_user1
    if user1(i)==0
        user1(i)=-1;
    end
end
for i=1:length_user2
    if user2(i)==0
        user2(i)=-1;
    end
end
fc=1; %%carrier frequency
eb=2; %% energy per bit
tb=1; %% time per bit of message sequence

%%% CDMA transmitter for user1
t=0.01:0.01:tb*length_user1; %0.01

%%plotting base band signal for user1
basebandsig1=[];
for i=1:length_user1
    for j=0.01:0.01:tb%0.01
          if user1(i)==1 
              basebandsig1=[basebandsig1 1];
          else 
              basebandsig1=[basebandsig1 -1];
          end
    end
end
% figure
% plot(basebandsig1)
% axis([0 100*length_user1 -1.5 1.5]);
% title('original binary sequence for user1 is')

%% BPSK MODULATION FOR USER 1
bpskmod1=[];
for i=1:length_user1
    for j=0.01:0.01:tb
     bpskmod1=[bpskmod1 sqrt(2*eb)*user1(i)*cos(2*pi*fc*j)];
 end
end
length(bpskmod1)

%  figure
%  plot(bpskmod1)
%  axis([0 100*length_user1 -2 2]);
%  title(' BPSK signal for user 1 is')
%% plot fft of BPSK sequence
% figure
% plot(real(fft(bpskmod1)))
% title('FFT of BPSK signal for user1 is')
%% PN generator for user1
%% let initial seed for user1 is 1000
seed1=[1 -1 1 -1];  %convert it into bipolar NRZ format 
spreadspectrum1=[];
pn1=[];
for i=1:length_user1
    for j=1:10 %chip rate is 10 times the bit rate
        pn1=[pn1 seed1(4)];  
        if seed1 (4)==seed1(3) temp=-1;
          else temp=1;
          end
              seed1(4)=seed1(3);
              seed1(3)=seed1(2);
              seed1(2)=seed1(1);
              seed1(1)=temp;
    end
         spreadspectrum1=[spreadspectrum1 user1(i)*pn1];
end


pnupsampled1=[];
len_pn1=length(pn1);
for i=1:len_pn1
    for j=0.1:0.1:tb
          if pn1(i)==1 
              pnupsampled1=[pnupsampled1 1];
          else 
              pnupsampled1=[pnupsampled1 -1];
          end
    end
end
length_pnupsampled1=length(pnupsampled1);
 sigtx1=bpskmod1.*pnupsampled1;
 
% figure
% plot(sigtx1)
% axis([0 100*length_user1 -2 2]);
% title(' spread spectrum signal transmitted for user 1 is')
% %% plot fft of BPSK sequence
% figure
% plot(real(fft(sigtx1)))
% title('FFT of spread spectrum signal transmitted for user1 is')
rxcode1=pnupsampled1.*sigtx1;
% figure
% plot(rxcode1)
% axis([0 100*length_user1 -2 2]);
% title(' spread spectrum signal transmitted for user 1 is')
%%% CDMA transmitter for user2
t=0.01:0.01:tb*length_user2; %0.01
%%plotting base band signal for user2
basebandsig2=[];
for i=1:length_user2
    for j=0.01:0.01:tb%0.01
          if user2(i)==1 
              basebandsig2=[basebandsig2 1];
          else 
              basebandsig2=[basebandsig2 -1];
          end
    end
end
% figure
% plot(basebandsig2)
% axis([0 100*length_user2 -1.5 1.5]);
% title('original binary sequence for user2 is')
%%%% BPSK MODULATION FOR USER 2
bpskmod2=[];
for i=1:length_user2
    for j=0.01:0.01:tb
     bpskmod2=[bpskmod2 sqrt(2*eb)*user2(i)*cos(2*pi*fc*j)];
 end
end
%  figure
%  plot(bpskmod2)
% axis([0 100*length_user2 -2 2]);
%  title(' BPSK signal for user 2 is')
%% plot fft of BPSK sequence
% figure
% plot(real(fft(bpskmod2)))
% title('FFT of BPSK signal for user2 is')
%% PN generator for user2
%% let initial seed for user2 is 1000
seed2=[-1 1 -1 1];  %convert it into bipolar NRZ format 
spreadspectrum2=[];
pn2=[];
for i=1:length_user2
    for j=1:10 %chip rate is 10 times the bit rate
        pn2=[pn2 seed2(4)];  
        if seed2 (4)==seed2(3) temp=-1;
          else temp=1;
          end
              seed2(4)=seed2(3);
              seed2(3)=seed2(2);
              seed2(2)=seed2(1);
              seed2(1)=temp;
    end
         spreadspectrum2=[spreadspectrum2 user2(i)*pn2];
end
pnupsampled2=[];
len_pn2=length(pn2);
for i=1:len_pn2
    for j=0.1:0.1:tb
          if pn2(i)==1 
              pnupsampled2=[pnupsampled2 1];
          else 
              pnupsampled2=[pnupsampled2 -1];
          end
    end
end
length_pnupsampled2=length(pnupsampled2);
 sigtx2=bpskmod2.*pnupsampled2;
 
% figure
% plot(sigtx2)
% axis([0 100*length_user2 -2 2]);
% title(' spread spectrum signal transmitted for user 2 is')
% %% plot fft of BPSK sequence
% figure
% plot(real(fft(sigtx2)))
% title('FFT of spread spectrum signal transmitted for user2 is')
 rxcode2=pnupsampled2.*sigtx2;
% figure
% plot(rxcode2)
% axis([0 100*length_user2 -2 2]);
% title(' spread spectrum signal transmitted for user 2 is')
%%signal by adding the above 2 signals%%
composite_signal=sigtx1+sigtx2;
%%AWGN CHANNEL%%
composite_signal=awgn(composite_signal,snr_in_dbs);  %% SNR of % dbs

%%DMODULATION FOR USER 1%%
rx1=composite_signal.*pnupsampled1;
% figure
% plot(rx1)
%% BPSK DEMODULATION FOR USER 1
demodcar1=[];
for i=1:length_user1
    for j=0.01:0.01:tb
     demodcar1=[demodcar1 sqrt(2*eb)*cos(2*pi*fc*j)];
 end
end
bpskdemod1=rx1.*demodcar1;
% figure
% plot(bpskdemod1)
% title('o/p of bpsk demod for user 1 is ')
len_dmod1=length(bpskdemod1);
sum=zeros(1,len_dmod1/100);
for i=1:len_dmod1/100
    for j=(i-1)*100+1:i*100
        sum(i)=sum(i)+bpskdemod1(j);
    end
end
sum;
  
 rxbits1=[];
 for i=1:length_user1
    if sum(i)>0
        rxbits1=[rxbits1 1];
    else
        rxbits1=[rxbits1 0];
    end
 end
rxbits1;

%1st Image reconstruction 
rxbits11=rxbits1';
output_1=uint8(rxbits11);
output_1=output_1(1:len1);
b1=reshape(output_1,8,N)';
dec_1=bi2de(b1,'left-msb');
im_output1=reshape(dec_1(1:N),size(in11,1),size(in11,2),size(in11,3));
im_output11=imresize(im_output1,10);
figure();
imshow(im_output11);title('Output1');

%%DMODULATION FOR USER 2%%
rx2=composite_signal.*pnupsampled2;
% figure
% plot(rx2)
%% BPSK DEMODULATION FOR USER 2
demodcar2=[];
for i=1:length_user2
    for j=0.01:0.01:tb
     demodcar2=[demodcar2 sqrt(2*eb)*cos(2*pi*fc*j)];
 end
end
bpskdemod2=rx2.*demodcar2
% figure
% plot(bpskdemod2)
% title('o/p of bpsk demod for user 2 is ')
len_dmod2=length(bpskdemod2);
sum=zeros(1,len_dmod1/100);
for i=1:len_dmod2/100
    for j=(i-1)*100+1:i*100
        sum(i)=sum(i)+bpskdemod2(j);
    end
end
sum;
  
 rxbits2=[];
 for i=1:length_user2
    if sum(i)>0
        rxbits2=[rxbits2 1];
    else
        rxbits2=[rxbits2 0];
    end
 end
rxbits2;

%2nd Image Reconstruction 
rxbits22=rxbits2';
output_2=uint8(rxbits22);
output_2=output_2(1:Z);
b2=reshape(output_2,8,M)';
dec_2=bi2de(b2,'left-msb');
im_output2=reshape(dec_2(1:M),size(in22,1),size(in22,2),size(in22,3));
figure();
im_output22=imresize(im_output2,10);
imshow(im_output22);title('Output2');

end

