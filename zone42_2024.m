function [beta,lamb,kap,k,nM,M,nphi,phi,nstiff,stiff,netf,kap_com,eps_bot,eps_top,eps_st_bot,eps_st_top,Matrix_sum_com,Matrix_sum_ten,Rebar_com,Rebar_ten, Efficiency] = zone42_2024(x,Mcr,epsilon_cr,phicr,EIcr,rho,n,zeta,xi,alpha,eta,kappa,omega,tau,mu)
beta = x/10;
k = -(-2*tau-2*mu*beta+2*mu*tau-eta*tau^2-rho*n*kappa*beta+1+2*eta*tau-eta-zeta*rho*n*beta^2+sqrt(-4*eta*tau*zeta*rho*n*beta^2*alpha-4*mu*tau*zeta*rho*n*beta^2*alpha+2*eta*tau^2*zeta*rho*n*beta^2*alpha-xi*beta^2+xi*beta^2*eta+2*rho^2*n^2*kappa*beta^3*zeta+2*xi*beta^4*zeta*rho*n+2*xi*beta^3*rho*n*kappa+rho^2*n^2*kappa^2*beta^2+zeta^2*rho^2*n^2*beta^4+xi*beta^2*eta*tau^2-2*xi*beta^2*mu*tau-2*xi*beta^2*eta*tau+4*tau*zeta*rho*n*beta^2*alpha+4*mu*beta^3*zeta*rho*n*alpha+2*eta*zeta*rho*n*beta^2*alpha-2*xi*beta^4*zeta*rho*n*alpha+2*xi*beta^3*mu+2*xi*beta^2*tau-2*zeta*rho*n*beta^2*alpha))/(2*tau-2*eta*tau-2*mu*tau+2*mu*beta+eta*tau^2+eta-1-xi*beta^2);
lamb = k*beta/(1-k);
kap =(-alpha+k)*beta/(k-1);
nM =((-2*beta^3*xi+2*eta*tau^3+3*beta^2*mu+(-3*mu-3*eta+3)*tau^2+eta-1)*k^3+(-6*zeta*rho*n*beta^3+(-6*kappa*n*rho-9*mu)*beta^2-6*eta*tau^3+(9*mu+9*eta-9)*tau^2-3*eta+3)*k^2+(-12*zeta*rho*n*(alpha-1)*beta^3+(6*kappa*n*(alpha+1)*rho+9*mu)*beta^2+6*eta*tau^3+(-9*mu-9*eta+9)*tau^2+3*eta-3)*k-6*zeta*n*rho*(alpha-1)^2*beta^3+(-6*alpha*kappa*n*rho-3*mu)*beta^2-2*eta*tau^3+(3*mu+3*eta-3)*tau^2-eta+1)/(beta^2*(k-1));
M = nM*Mcr;
nphi = (1/2)*beta/(1-k);
phi = nphi*phicr;
nstiff =-(4*(-beta^3*xi+(3/2*(beta^2))*mu+eta*tau^3+(-3/2*mu-3/2*eta+3/2)*tau^2-1/2+(1/2)*eta))*k^3/beta^3-(4*(-3*zeta*rho*n*beta^3+(-9/2*mu-3*rho*n*kappa)*beta^2-3*eta*tau^3+(9/2*mu+9/2*eta-9/2)*tau^2+3/2-3/2*eta))*k^2/beta^3-(4*(-6*zeta*rho*n*(alpha-1)*beta^3+(3*kappa*n*(alpha+1)*rho+9/2*mu)*beta^2+3*eta*tau^3+(-9/2*mu-9/2*eta+9/2)*tau^2-3/2+3/2*eta))*k/beta^3-(4*(-3*zeta*n*rho*(alpha-1)^2*beta^3+(-3/2*mu-3*alpha*kappa*n*rho)*beta^2-eta*tau^3+(-3/2+3/2*mu+3/2*eta)*tau^2+1/2-(1/2)*eta))/beta^3;
stiff = nstiff*EIcr;
netf = (1/2)*((eta*tau^2+(-2*mu-2*eta+2)*tau-beta^2*xi+2*beta*mu+eta-1)*k^2+(-2*eta*tau^2+(4*mu+4*eta-4)*tau-2*zeta*rho*n*beta^2+(-2*kappa*n*rho-4*mu)*beta-2*eta+2)*k+eta*tau^2+(-2*mu-2*eta+2)*tau-2*zeta*rho*n*(alpha-1)*beta^2+(2*kappa*n*rho+2*mu)*beta+eta-1)/((k-1)*beta);
kap_com = -(k-1+alpha)*beta/(k-1);
eps_bot = beta*epsilon_cr;
eps_top = lamb*epsilon_cr;
eps_st_bot = kap*epsilon_cr;
eps_st_top = kap_com*epsilon_cr;
Matrix_sum_com =-k^2*beta*epsilon_cr/(2*k - 2);
Matrix_sum_ten =(k - 1)*(-eta*tau^2/2 + (mu + eta - 1)*tau - mu*beta - eta/2 + 1/2)*epsilon_cr/beta;
Rebar_com =-zeta*rho*n*(k - 1 + alpha)*beta*epsilon_cr/(k - 1);
Rebar_ten =rho*n*kappa*epsilon_cr;

Efficiency = [k^2*xi/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), 2*n*rho*(k - 1 + alpha)*zeta/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), 2*(k - 1)*(-eta*tau^2/2 + (mu + eta - 1)*tau - mu*beta - eta/2 + 1/2)/(-eta*(k - 1)*tau^2 + 2*(mu + eta - 1)*(k - 1)*tau + (1 - k)*eta + (-2*beta*mu + 1)*k + 2*beta*kappa*n*rho + 2*mu*beta - 1), 2*n*rho*beta*kappa/(-eta*(k - 1)*tau^2 + 2*(mu + eta - 1)*(k - 1)*tau + (-2*beta*mu - eta + 1)*k + 2*beta*kappa*n*rho + 2*mu*beta + eta - 1)];
end