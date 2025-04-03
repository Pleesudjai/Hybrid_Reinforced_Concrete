function [beta,lamb,kap,k,nM,M,nphi,phi,nstiff,stiff,netf,kap_com,eps_bot,eps_top,eps_st_bot,eps_st_top,Matrix_sum_com,Matrix_sum_ten,Rebar_com,Rebar_ten,Efficiency] = zone31_2024(x,Mcr,epsilon_cr,phicr,EIcr,rho,n,zeta,xi,alpha,eta,kappa,omega)
beta = x/10;
k =-(-zeta*rho*n*beta^2-rho*n*beta^2-beta^2*eta-beta*omega*xi-xi*omega^2+2*eta*beta+sqrt(beta^4*n^2*rho^2*zeta^2+2*alpha*beta^4*eta*n*rho*zeta+4*alpha*beta^3*n*omega*rho*xi*zeta+2*alpha*beta^2*n*omega^2*rho*xi*zeta+2*beta^4*n^2*rho^2*zeta-2*alpha*beta^4*eta*n*rho-4*alpha*beta^3*eta*n*rho*zeta-4*alpha*beta^3*n*omega*rho*xi-2*alpha*beta^2*n*omega^2*rho*xi+beta^4*n^2*rho^2-2*beta^3*n*omega*rho*xi*zeta+4*alpha*beta^3*eta*n*rho+4*alpha*beta^3*n*rho*zeta+2*alpha*beta^2*eta*n*rho*zeta+2*beta^4*eta*n*rho+2*beta^3*n*omega*rho*xi+2*beta^2*n*omega^2*rho*xi-4*alpha*beta^3*n*rho-2*alpha*beta^2*eta*n*rho-2*alpha*beta^2*n*rho*zeta-4*beta^3*eta*n*rho+beta^2*omega^2*xi^2+2*alpha*beta^2*n*rho+4*beta^3*n*rho+2*beta^2*eta*n*rho-2*beta^2*n*rho)-2*beta-eta+1)/(beta^2*eta+2*beta*omega*xi+omega^2*xi-2*beta*eta+2*beta+eta-1) ;
lamb = k*beta/(1-k);
kap =(-alpha+k)*beta/(k-1);
nM =((2*beta^3*eta+(3*omega*xi-3*eta+3)*beta^2-xi*omega^3+eta-1)*k^3+((-6*eta+(-6*zeta-6)*n*rho)*beta^3+(-3*omega*xi+9*eta-9)*beta^2+3*xi*omega^3-3*eta+3)*k^2+((6*eta+((-12*alpha+12)*zeta+12*alpha)*n*rho)*beta^3+(-9*eta+9)*beta^2-3*xi*omega^3+3*eta-3)*k+(-2*eta+((-6*alpha^2+12*alpha-6)*zeta-6*alpha^2)*n*rho)*beta^3+(3*eta-3)*beta^2+xi*omega^3-eta+1)/((k-1)*beta^2) ;
M = nM*Mcr;
nphi = (1/2)*beta/(1-k);
phi = nphi*phicr;
nstiff =(12*(-(1/3)*beta^3*eta+(-1/2+(1/2)*eta-(1/2)*xi*omega)*beta^2-(1/6)*eta+1/6+(1/6)*xi*omega^3))*k^3/beta^3+(12*((eta+(zeta+1)*n*rho)*beta^3+(3/2-3/2*eta+(1/2)*xi*omega)*beta^2+(1/2)*eta-(1/2)*xi*omega^3-1/2))*k^2/beta^3+(12*((-eta+((2*alpha-2)*zeta-2*alpha)*n*rho)*beta^3+(-3/2+3/2*eta)*beta^2+1/2-(1/2)*eta+(1/2)*xi*omega^3))*k/beta^3+(12*(((1/3)*eta+((alpha^2-2*alpha+1)*zeta+alpha^2)*n*rho)*beta^3+(-(1/2)*eta+1/2)*beta^2-(1/6)*xi*omega^3-1/6+(1/6)*eta))/beta^3;
stiff = nstiff*EIcr;
netf =(1/2)*((beta^2*eta+(2*(omega*xi-eta+1))*beta+xi*omega^2+eta-1)*k^2+((-2*eta+(-2*zeta-2)*n*rho)*beta^2+(2*(-omega*xi+2*eta-2))*beta-2*xi*omega^2-2*eta+2)*k+(eta+((-2*alpha+2)*zeta+2*alpha)*n*rho)*beta^2+(2*(-eta+1))*beta+xi*omega^2+eta-1)/((k-1)*beta) ;
kap_com = -(k-1+alpha)*beta/(k-1);
eps_bot = beta*epsilon_cr;
eps_top = lamb*epsilon_cr;
eps_st_bot = kap*epsilon_cr;
eps_st_top = kap_com*epsilon_cr;
Matrix_sum_com =(epsilon_cr*(k - 1)*omega^2 + 2*epsilon_cr*k*beta*omega)/(2*beta);
Matrix_sum_ten =-((k - 1)*epsilon_cr*(beta^2*eta + (-2*eta + 2)*beta + eta - 1))/(2*beta);
Rebar_com =-zeta*rho*n*(k - 1 + alpha)*beta*epsilon_cr/(k - 1);
Rebar_ten =rho*n*(-alpha + k)*beta*epsilon_cr/(k - 1);

Efficiency = [(xi*omega^2*(k - 1)^2 + 2*omega*k*xi*(k - 1)*beta - 2*zeta*rho*n*(k - 1 + alpha)*beta^2)/(xi*omega^2*(k - 1)^2 + 2*omega*k*xi*(k - 1)*beta), -2*n*rho*beta^2*zeta*(k - 1 + alpha)/(xi*omega^2*(k - 1)^2 + 2*omega*k*xi*(k - 1)*beta), ((eta*(k - 1)^2 - 2*n*rho*(-alpha + k))*beta^2 - 2*(k - 1)^2*(eta - 1)*beta + (k - 1)^2*(eta - 1))/((k - 1)^2*(beta^2*eta + (-2*eta + 2)*beta + eta - 1)), -2*rho*n*(-alpha + k)*beta^2/((k - 1)^2*(beta^2*eta + (-2*eta + 2)*beta + eta - 1))];

end