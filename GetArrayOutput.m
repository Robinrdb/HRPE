function X = GetArrayOutput(M,mu,SNR,N,rho)
% X = GetArrayOutput(M,mu,SNR,N[,rho])
% Computes the output matrix of a uniform linear array.
% Uncorrelated unity power Gaussian signals are used
% the noise is additive white and Gaussian.
%
% Input:
%   M   : Number of antenna elements
%   mu  : Vector of spatial frequencies
%   SNR : Signal-to-noise ratio in dB
%   N   : Number of snapshots
%   rho : the correlation coefficient (optional, defaults to 0)
%
% Output:
%   X   : Array output matrix. Each column contains a snapshot

    if (nargin<5), rho = 0; end;
    if (nargin<4), error('GetArrayOutput::Too few parameters'); end;
    mu = mu(:);
    d  = length(mu);            % number of signals
    A  = exp(j*(0:M-1)'*mu');   % array steering matrix
    S  = (randn(d,N) + j * randn(d,N)) / sqrt(2);  % signal matrix
    Rs = ones(d)*rho + eye(d)*(1-rho);
    S  = sqrtm(Rs)*S;
    NoisePwr = 10^(-SNR/10);
    N  = sqrt(NoisePwr) * (randn(M,N) + j * randn(M,N)) / sqrt(2); % noise 
    X  = A * S + N;             % array output
