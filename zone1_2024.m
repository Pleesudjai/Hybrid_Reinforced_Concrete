function [beta,lamb,kap,k,nM,M,nphi,phi,nstiff,stiff,netf,kap_com,eps_bot,eps_top,eps_st_bot,eps_st_top,Matrix_sum_com,Matrix_sum_ten,Rebar_com,Rebar_ten,Efficiency] = zone1_2024(x,Mcr,epsilon_cr,phicr,EIcr,rho,n,zeta,xi,alpha)
beta = x/10;
k =(sqrt(n^2*(zeta+1)^2*rho^2-2*(((xi-1)*alpha-xi)*zeta-1+(-xi+1)*alpha)*n*rho+xi)-1+(-zeta-1)*n*rho)/(xi-1);
lamb = k*beta/(1-k);
kap =(-alpha+k)*beta/(k-1);
nM =-2*beta*((xi-1)*k^3+(3+3*n*(zeta+1)*rho)*k^2+(-3+6*n*((alpha-1)*zeta-alpha)*rho)*k+1+3*n*((alpha-1)^2*zeta+alpha^2)*rho)/(k-1);
M = nM*Mcr;
nphi = (1/2)*beta/(1-k);
phi = nphi*phicr;
nstiff =(4*xi-4)*k^3+(12*n*(zeta+1)*rho+12)*k^2+((12*((2*alpha-2)*zeta-2*alpha))*n*rho-12)*k+12*n*((alpha-1)^2*zeta+alpha^2)*rho+4;
stiff = nstiff*EIcr;
netf = -beta*((xi-1)*k^2+(2+2*n*(zeta+1)*rho)*k-1+2*n*((alpha-1)*zeta-alpha)*rho)/(2*k-2);
kap_com = -(k-1+alpha)*beta/(k-1);
eps_bot = beta*epsilon_cr;
eps_top = lamb*epsilon_cr;
eps_st_bot = kap*epsilon_cr;
eps_st_top = kap_com*epsilon_cr;
Matrix_sum_com =-k^2*beta*epsilon_cr/(2*k - 2);
Matrix_sum_ten =-((k - 1)*epsilon_cr*beta)/2;
Rebar_com =-zeta*rho*n*(k - 1 + alpha)*beta*epsilon_cr/(k - 1);
Rebar_ten =rho*n*(-alpha + k)*beta*epsilon_cr/(k - 1);

Efficiency = [k^2*xi/(2*zeta*rho*n*(k - 1 + alpha) + xi*k^2), 2*n*rho*(k - 1 + alpha)*zeta/(2*zeta*rho*n*(k - 1 + alpha) + xi*k^2), (k - 1)^2/(k^2 + (-2*n*rho - 2)*k + 2*rho*n*alpha + 1), -2*(-alpha + k)*n*rho/(k^2 + (-2*n*rho - 2)*k + 2*rho*n*alpha + 1)];
end