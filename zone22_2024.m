function [beta,lamb,kap,k,nM,M,nphi,phi,nstiff,stiff,netf,kap_com,eps_bot,eps_top,eps_st_bot,eps_st_top,Matrix_sum_com,Matrix_sum_ten,Rebar_com,Rebar_ten, Efficiency] = zone22_2024(x,Mcr,epsilon_cr,phicr,EIcr,rho,n,zeta,xi,alpha,eta,kappa)
beta = x/10;
k =-(-zeta*rho*n*beta^2-rho*n*kappa*beta-beta^2*eta+2*eta*beta+sqrt(beta^4*n^2*rho^2*zeta^2+2*alpha*beta^4*eta*n*rho*zeta-2*alpha*beta^4*n*rho*xi*zeta+2*beta^3*kappa*n^2*rho^2*zeta-4*alpha*beta^3*eta*n*rho*zeta+2*beta^4*n*rho*xi*zeta+beta^2*kappa^2*n^2*rho^2+4*alpha*beta^3*n*rho*zeta+2*alpha*beta^2*eta*n*rho*zeta+2*beta^3*kappa*n*rho*xi-2*alpha*beta^2*n*rho*zeta+beta^4*eta*xi-2*beta^3*eta*xi+2*beta^3*xi+beta^2*eta*xi-beta^2*xi)-2*beta-eta+1)/(beta^2*eta-beta^2*xi-2*beta*eta+2*beta+eta-1);
lamb = k*beta/(1-k);
kap =(-alpha+k)*beta/(k-1);
nM = (((-2*xi+2*eta)*beta^3-(3*(eta-1))*beta^2+eta-1)*k^3+((-6*n*rho*zeta-6*eta)*beta^3-(3*(2*kappa*n*rho-3*eta+3))*beta^2-3*eta+3)*k^2+((6*eta+(-12*alpha+12)*zeta*n*rho)*beta^3-(3*(-2*alpha*kappa*n*rho-2*kappa*n*rho+3*eta-3))*beta^2+3*eta-3)*k+(-2*eta+(-6*alpha^2+12*alpha-6)*zeta*n*rho)*beta^3-(3*(2*alpha*kappa*n*rho-eta+1))*beta^2-eta+1)/((k-1)*beta^2);
M = nM*Mcr;
nphi = (1/2)*beta/(1-k);
phi = nphi*phicr;
nstiff =(12*(((1/3)*xi-(1/3)*eta)*beta^3-(-(1/2)*eta+1/2)*beta^2-(1/6)*eta+1/6))*k^3/beta^3+(12*((n*rho*zeta+eta)*beta^3-(3/2*eta-3/2-n*rho*kappa)*beta^2+(1/2)*eta-1/2))*k^2/beta^3+(12*((-eta+(2*alpha-2)*zeta*n*rho)*beta^3-(n*rho*kappa-3/2*eta+3/2+alpha*kappa*n*rho)*beta^2-(1/2)*eta+1/2))*k/beta^3+(12*(((1/3)*eta+(alpha^2-2*alpha+1)*zeta*n*rho)*beta^3-(-alpha*kappa*n*rho+(1/2)*eta-1/2)*beta^2+(1/6)*eta-1/6))/beta^3;
stiff = nstiff*EIcr;
netf = (1/2)*(((-xi+eta)*beta^2-(2*(eta-1))*beta+eta-1)*k^2+((-2*n*rho*zeta-2*eta)*beta^2-(2*(kappa*n*rho-2*eta+2))*beta-2*eta+2)*k+(eta+(-2*alpha+2)*zeta*n*rho)*beta^2-(2*(-kappa*n*rho+eta-1))*beta+eta-1)/((k-1)*beta);
kap_com = -(k-1+alpha)*beta/(k-1);
eps_bot = beta*epsilon_cr;
eps_top = lamb*epsilon_cr;
eps_st_bot = kap*epsilon_cr;
eps_st_top = kap_com*epsilon_cr;
Matrix_sum_com =-k^2*beta*epsilon_cr/(2*k - 2);
Matrix_sum_ten =-((k - 1)*epsilon_cr*(beta^2*eta + (-2*eta + 2)*beta + eta - 1))/(2*beta);
Rebar_com =-zeta*rho*n*(k - 1 + alpha)*beta*epsilon_cr/(k - 1);
Rebar_ten =rho*n*kappa*epsilon_cr;

Efficiency = [(k^2*xi/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi)), (2*n*rho*(k - 1 + alpha)*zeta/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi)), ((beta^2*eta + (-2*eta + 2)*beta + eta - 1)*(k - 1)/(eta*(k - 1)*beta^2 + ((-2*k + 2)*eta - 2*n*rho*kappa + 2*k - 2)*beta + (eta - 1)*(k - 1))), (-2*rho*n*kappa*beta/(eta*(k - 1)*beta^2 + ((-2*k + 2)*eta - 2*n*rho*kappa + 2*k - 2)*beta + (eta - 1)*(k - 1)))];

end