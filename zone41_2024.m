function [beta,lamb,kap,k,nM,M,nphi,phi,nstiff,stiff,netf,kap_com,eps_bot,eps_top,eps_st_bot,eps_st_top,Matrix_sum_com,Matrix_sum_ten,Rebar_com,Rebar_ten,Efficiency] = zone41_2024(x,Mcr,epsilon_cr,phicr,EIcr,rho,n,zeta,xi,alpha,eta,kappa,omega,tau,mu)
beta = x/10;
k =(-(-1+tau)^2*eta-n*rho*(zeta+1)*beta^2-2*beta*mu+(2*mu-2)*tau+sqrt(beta^2*(n^2*beta^2*(zeta+1)^2*rho^2-4*n*((-(1/2)*eta*tau^2+(mu+eta-1)*tau+(1/2)*beta^2*xi-beta*mu-(1/2)*eta+1/2)*(zeta-1)*alpha-(1/2)*eta*tau^2+(mu+eta-1)*tau-(1/2)*beta^2*xi*zeta-beta*mu-(1/2)*eta+1/2)*rho-2*xi*(-(1/2)*eta*tau^2+(mu+eta-1)*tau-beta*mu-(1/2)*eta+1/2)))+1)/(-(-1+tau)^2*eta+beta^2*xi-2*beta*mu+1+(2*mu-2)*tau);
lamb = k*beta/(1-k);
kap =(-alpha+k)*beta/(k-1);
nM =((2*eta*tau^3+(-3*mu-3*eta+3)*tau^2-2*beta^3*xi+3*beta^2*mu+eta-1)*k^3+(-6*eta*tau^3+(9*mu+9*eta-9)*tau^2-6*n*rho*(zeta+1)*beta^3-9*beta^2*mu-3*eta+3)*k^2+(6*eta*tau^3+(-9*mu-9*eta+9)*tau^2-12*n*((alpha-1)*zeta-alpha)*rho*beta^3+9*beta^2*mu+3*eta-3)*k-2*eta*tau^3+(3*mu+3*eta-3)*tau^2-(6*((alpha-1)^2*zeta+alpha^2))*n*rho*beta^3-3*beta^2*mu-eta+1)/(beta^2*(k-1));
M = nM*Mcr;
nphi = (1/2)*beta/(1-k);
phi = nphi*phicr;
nstiff =-(4*(eta*tau^3+(-3/2*mu-3/2*eta+3/2)*tau^2-1/2-beta^3*xi+(3/2*(beta^2))*mu+(1/2)*eta))*k^3/beta^3-(4*(-3*eta*tau^3+(9/2*mu+9/2*eta-9/2)*tau^2-3*n*rho*(zeta+1)*beta^3-(9/2*(beta^2))*mu+3/2-3/2*eta))*k^2/beta^3-(4*(3*eta*tau^3+(-9/2*mu-9/2*eta+9/2)*tau^2-6*n*((alpha-1)*zeta-alpha)*rho*beta^3+(9/2*(beta^2))*mu-3/2+3/2*eta))*k/beta^3-(4*(-eta*tau^3+(-3/2+3/2*mu+3/2*eta)*tau^2-(3*((alpha-1)^2*zeta+alpha^2))*n*rho*beta^3-(3/2*(beta^2))*mu+1/2-(1/2)*eta))/beta^3;
stiff = nstiff*EIcr;
netf = (1/2)*((eta*tau^2+(-2*mu-2*eta+2)*tau-beta^2*xi+2*beta*mu+eta-1)*k^2+(-2*eta*tau^2+(4*mu+4*eta-4)*tau-2*n*rho*(zeta+1)*beta^2-4*beta*mu-2*eta+2)*k+eta*tau^2+(-2*mu-2*eta+2)*tau-2*n*((alpha-1)*zeta-alpha)*rho*beta^2+2*beta*mu+eta-1)/((k-1)*beta);
kap_com = -(k-1+alpha)*beta/(k-1);
eps_bot = beta*epsilon_cr;
eps_top = lamb*epsilon_cr;
eps_st_bot = kap*epsilon_cr;
eps_st_top = kap_com*epsilon_cr;
Matrix_sum_com =-k^2*beta*epsilon_cr/(2*k - 2);
Matrix_sum_ten =(k - 1)*(-eta*tau^2/2 + (mu + eta - 1)*tau - mu*beta - eta/2 + 1/2)*epsilon_cr/beta;
Rebar_com =-zeta*rho*n*(k - 1 + alpha)*beta*epsilon_cr/(k - 1);
Rebar_ten =rho*n*(-alpha + k)*beta*epsilon_cr/(k - 1);

Efficiency = [k^2*xi/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), 2*n*rho*(k - 1 + alpha)*zeta/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), 2*(-eta*tau^2/2 + (mu + eta - 1)*tau - mu*beta - eta/2 + 1/2)*(k - 1)^2/(-eta*(k - 1)^2*tau^2 + 2*(k - 1)^2*(mu + eta - 1)*tau + (-2*beta*mu - eta + 1)*k^2 + (2*beta^2*n*rho + 4*beta*mu + 2*eta - 2)*k - 2*rho*n*beta^2*alpha - 2*mu*beta - eta + 1), 2*rho*n*(-alpha + k)*beta^2/((-eta*tau^2 + (2*mu + 2*eta - 2)*tau - 2*mu*beta - eta + 1)*k^2 + (2*eta*tau^2 + (-4*mu - 4*eta + 4)*tau + 2*n*rho*beta^2 + 4*mu*beta + 2*eta - 2)*k - eta*tau^2 + (2*mu + 2*eta - 2)*tau - 2*rho*n*beta^2*alpha - 2*mu*beta - eta + 1)];
end