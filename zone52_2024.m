function [beta,lamb,kap,k,nM,M,nphi,phi,nstiff,stiff,netf,kap_com,eps_bot,eps_top,eps_st_bot,eps_st_top,Matrix_sum_com,Matrix_sum_ten,Rebar_com,Rebar_ten,Efficiency ] = zone52_2024(x,Mcr,epsilon_cr,phicr,EIcr,rho,n,zeta,xi,alpha,eta,kappa,omega,tau,mu)
beta = x/10;
k =-(-zeta*rho*n*beta^2-rho*n*kappa*beta-beta*omega*xi-eta*tau^2-xi*omega^2-2*beta*mu+2*eta*tau+2*mu*tau+sqrt(beta^4*n^2*rho^2*zeta^2+4*alpha*beta^3*n*omega*rho*xi*zeta+2*alpha*beta^2*eta*n*rho*tau^2*zeta+2*alpha*beta^2*n*omega^2*rho*xi*zeta+2*beta^3*kappa*n^2*rho^2*zeta+4*alpha*beta^3*mu*n*rho*zeta-4*alpha*beta^2*eta*n*rho*tau*zeta-4*alpha*beta^2*mu*n*rho*tau*zeta-2*beta^3*n*omega*rho*xi*zeta+beta^2*kappa^2*n^2*rho^2+2*alpha*beta^2*eta*n*rho*zeta+4*alpha*beta^2*n*rho*tau*zeta-2*beta^2*kappa*n*omega*rho*xi-2*alpha*beta^2*n*rho*zeta+beta^2*omega^2*xi^2)-eta-2*tau+1)/(2*beta*omega*xi+eta*tau^2+omega^2*xi+2*beta*mu-2*eta*tau-2*mu*tau+eta+2*tau-1);
lamb = k*beta/(1-k);
kap =(-alpha+k)*beta/(k-1);
nM = (((3*omega*xi+3*mu)*beta^2+2*eta*tau^3+(-3*mu-3*eta+3)*tau^2-xi*omega^3+eta-1)*k^3+(-6*zeta*rho*n*beta^3+(-6*kappa*n*rho-3*omega*xi-9*mu)*beta^2-6*eta*tau^3+(9*mu+9*eta-9)*tau^2+3*xi*omega^3-3*eta+3)*k^2+(-12*zeta*n*rho*(alpha-1)*beta^3+(6*kappa*n*(alpha+1)*rho+9*mu)*beta^2+6*eta*tau^3+(-9*mu-9*eta+9)*tau^2-3*xi*omega^3+3*eta-3)*k-6*zeta*n*rho*(alpha-1)^2*beta^3+(-6*alpha*kappa*n*rho-3*mu)*beta^2-2*eta*tau^3+(3*mu+3*eta-3)*tau^2+xi*omega^3-eta+1)/(beta^2*(-1+k));
M = nM*Mcr;
nphi = (1/2)*beta/(1-k);
phi = nphi*phicr;
nstiff =-(4*((3/2*mu+(3/2*xi)*omega)*beta^2+eta*tau^3+(3/2-3/2*mu-3/2*eta)*tau^2-1/2+(1/2)*eta-(1/2)*xi*omega^3))*k^3/beta^3-(4*(-3*zeta*rho*n*beta^3+(-3*kappa*n*rho-9/2*mu-(3/2*xi)*omega)*beta^2-3*eta*tau^3+(-9/2+9/2*mu+9/2*eta)*tau^2+(3/2*xi)*omega^3-3/2*eta+3/2))*k^2/beta^3-(4*(-6*zeta*n*rho*(alpha-1)*beta^3+(3*kappa*n*(alpha+1)*rho+9/2*mu)*beta^2+3*eta*tau^3+(9/2-9/2*mu-9/2*eta)*tau^2-3/2-(3/2*xi)*omega^3+3/2*eta))*k/beta^3-(4*(-3*zeta*n*rho*(alpha-1)^2*beta^3+(-3*alpha*kappa*n*rho-3/2*mu)*beta^2-eta*tau^3+(3/2*mu+3/2*eta-3/2)*tau^2+1/2+(1/2)*xi*omega^3-(1/2)*eta))/beta^3;
stiff = nstiff*EIcr;
netf= (1/2)*((eta*tau^2+(-2*mu-2*eta+2)*tau+(2*omega*xi+2*mu)*beta+xi*omega^2+eta-1)*k^2+(-2*eta*tau^2+(4*mu+4*eta-4)*tau-2*zeta*rho*n*beta^2+(-2*kappa*n*rho-2*omega*xi-4*mu)*beta-2*xi*omega^2-2*eta+2)*k+eta*tau^2+(-2*mu-2*eta+2)*tau-2*zeta*n*rho*(alpha-1)*beta^2+(2*kappa*n*rho+2*mu)*beta+xi*omega^2+eta-1)/((-1+k)*beta);
kap_com = -(k-1+alpha)*beta/(k-1);
eps_bot = beta*epsilon_cr;
eps_top = lamb*epsilon_cr;
eps_st_bot = kap*epsilon_cr;
eps_st_top = kap_com*epsilon_cr;
Matrix_sum_com =(epsilon_cr*(k - 1)*omega^2 + 2*epsilon_cr*k*beta*omega)/(2*beta);
Matrix_sum_ten =(k - 1)*(-eta*tau^2/2 + (mu + eta - 1)*tau - mu*beta - eta/2 + 1/2)*epsilon_cr/beta;
Rebar_com =-zeta*rho*n*(k - 1 + alpha)*beta*epsilon_cr/(k - 1);
Rebar_ten =rho*n*kappa*epsilon_cr;

Efficiency = [(xi*omega^2*(k - 1)^2 + 2*omega*k*xi*(k - 1)*beta)/(xi*omega^2*(k - 1)^2 + 2*omega*k*xi*(k - 1)*beta - 2*zeta*rho*n*(k - 1 + alpha)*beta^2), -2*n*rho*beta^2*zeta*(k - 1 + alpha)/(xi*omega^2*(k - 1)^2 + 2*omega*k*xi*(k - 1)*beta - 2*zeta*rho*n*(k - 1 + alpha)*beta^2), 2*(k - 1)*(-eta*tau^2/2 + (mu + eta - 1)*tau - mu*beta - eta/2 + 1/2)/(-eta*(k - 1)*tau^2 + 2*(mu + eta - 1)*(k - 1)*tau + (1 - k)*eta + (-2*beta*mu + 1)*k + 2*beta*kappa*n*rho + 2*mu*beta - 1), 2*n*rho*beta*kappa/(-eta*(k - 1)*tau^2 + 2*(mu + eta - 1)*(k - 1)*tau + (-2*beta*mu - eta + 1)*k + 2*beta*kappa*n*rho + 2*mu*beta + eta - 1)];
end