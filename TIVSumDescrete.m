
function [dy, sigma]= TIVSumDescrete(t,x,u,k1,k2,k3,kL,Ls,L0,M,b,alpha,sigma_i,G,varargin)


rho_m = x(1);    
rho_f = x(2);   


%% Evolution of Mean free path
%dL = -kL*(L-Ls);

L = Ls+(L0-Ls)*exp(-kL*u);
du = 

%% Evolution of mobile dislocation density
drho_m = (M/b)*(1/Ls - 1/L);

%% Evolution of forest density
drho_f = M*(k1*sqrt(rho_f))-k2*rho_f+k3*rho_m;   %M*((k1/(b*L))- k2*rhof);

%% output preparation
dy =  [drho_m;drho_f];
sigma = sigma_i + M*alpha*G*b*sqrt(rho_f);
end



