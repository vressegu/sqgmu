function t = fct_unity_approx5(N_t)
% XP must be a 2 x n matrix
% the result is a vector of size n
%

slop_size_ratio=6;
% slop_size_ratio=ceil(6*(N_t/1000));
% N_t=1000;

t=ones(1,N_t);
P_t=N_t/2;
sslop=ceil(N_t/slop_size_ratio);
t((P_t-sslop+1):P_t)= (-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;
t((P_t+2):(P_t+1+sslop))= (tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
% t(1:sslop)=(tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
% t(end-sslop+1:end)=(-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;

t(P_t+1)=0;

