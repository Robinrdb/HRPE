%%Intialization
%Snapshot number
N=500;
%element number
M=16;
%uncorrelated impinging signals number
d=4;
%directions of d impinging signals
mu_desired=[-0.6;-0.2;0.1;0.3];
%SNR in dB
SNR_matric=1:1:50;
%the correlation coefficient
rho=0;
errors1 = zeros(length(SNR_matric), NumberOfRuns);
errors2 = zeros(length(SNR_matric), NumberOfRuns);

%% ITERATION SNR
iter=0;
for i=1:50
    SNR=SNR_matric(i);
    for j=1:NumberOfRuns
%% GetArrayOutput
X = GetArrayOutput(M,mu_desired,SNR,N,rho);
%% Forward-Backward-Averaging
Z = [X,flip(eye(M))*conj(X)*flip(eye(N))];
%% Estimation of Signal Parameters via Rotational Invariance Techniques(ESPRIT)
mu1 = Esprit(X,d);
mu2 = Esprit(Z,d);
%% RMSE
errors1(i,j)=norm(mu_desired(:)-mu1(:))^2;
errors2(i,j)=norm(mu_desired(:)-mu2(:))^2;
    end
    iter=iter+1
end
semilogy(SNR_matric, mean(errors1,2))
hold on
semilogy(SNR_matric, mean(errors2,2))
grid on
legend('Without FBA','With FBA');
xlabel('SNR');ylabel('Error');
title('Estimation Error VS SNR in rho=0.99');
