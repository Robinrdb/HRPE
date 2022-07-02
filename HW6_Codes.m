clear;

%% Initialization

%ULA antenna numbers M
M=16;
% numbers of impinging signals
d=3;
% numbers of snapshots
N=500;
%standard beamwidth mub
mub=2*pi/M;
mu=mub.*[-3,0,3];
%SNR
SNR=10;
% channel parameters
S = ( sign(randn(d,N)) + j * sign(randn(d,N)) ) / sqrt(2);
W = ( randn(M,N) + j * randn(M,N) ) / sqrt(2) * 10^(-SNR/20);
A=zeros(M,d);
c=(0:M-1);
%a=amu(1,:) represents first row
for i=1:d
    A(:,i)=exp(j*c*mu(i));
end
%Input of filters
X=zeros(M,N);
X=A*S+W;
Rxx_est_new=0;
for i=1:N
Rxx_est_current=X(:,i)*X(:,i)';
Rxx_est_new=Rxx_est_current+Rxx_est_new;
end
Rxx_est=Rxx_est_new/N;
%% Sampling
NS=1000;
mus=-pi:2*pi/(NS-1):pi;
%% MVDR
S_MVDR=zeros(M,NS);
S_MVDR_mus=zeros(1,NS);
for i=1:NS
    S_MVDR(:,i)=exp(j*c*mus(i));
    S_MVDR_mus(i)=1/(S_MVDR(:,i)'*inv(Rxx_est)*(S_MVDR(:,i)));
end
hold on
semilogy(mus,normalize(S_MVDR_mus))
%% DFT
d=0;
f=0;
for n=1:N
    for k=0:M-1
        d1=X(k+1,n)*exp(-j*mus*k);
        d=d1+d;
    end
    f1=abs(d).*abs(d);
    f=f1+f;
end
S_DFT=f/(N*M);
hold on
semilogy(mus,normalize(S_DFT));
%% MUSIC
U0=null(A');
S_MUSIC=zeros(M,NS);
S_MUSIC_mus=zeros(1,NS);
for i=1:NS
     S_MUSIC(:,i)=exp(j*c*mus(i));
     S_MUSIC_mus(i)=(S_MUSIC(:,i)'*S_MUSIC(:,i))/(S_MUSIC(:,i)'*U0*U0'*(S_MUSIC(:,i)));
end
semilogy(mus,normalize(S_MUSIC_mus))
legend('MVDR','DFT','MUSIC');
xlabel('mu');ylabel('power spectra');
ylim([-2 10]);
title('Estimation of Direction of Desired Signal with (i)');
