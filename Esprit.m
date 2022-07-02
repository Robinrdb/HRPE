function mu = Esprit(X,d)
%Input:
%where X is the M × N matrix containing the array measurements 
%and d is the number of signal directions to look for
%
%Output:
%The return value mu is a vector containing the estimated
%spatial frequencies of the d wavefronts, sorted in a ascending order.
%% Selection Matrices
M=16;
N=500;
%subarray number small m(d+1<=m<M) when m=M-1 we have maximum overlap
m=M-1;
B=zeros(1,m).';
J1=[eye(m),B];
J2=[B,eye(m)];
[U,S,V] = svd(X);
Us=U(:,1:d);
%% Least Square Solution
pphi_LS=(pinv(J1*Us))*(J2*Us);
phi_LS1=eig(pphi_LS);
phi_LS=sort(phi_LS1);
mu=angle(phi_LS);
mu=sort(mu);