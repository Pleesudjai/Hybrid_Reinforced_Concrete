function [Output_Final,Output,casee] = Intersection_points_Tri_HRC(Mcr,phicr,rho,n,zeta,xi,alpha,eta,kappa,omega,tau,mu,epsilon_cr,d,E,b)

h=d;

%........................
% shared stage 0
%........................
Stage_0(1,1)= 0;
Stage_0(1,2)= 0;
Stage_0(1,3) =0 ;
Stage_0(1,4) =0 ;
Stage_0(1,5) =0 ;
Stage_0(1,6) =0 ;
beta_121 = 0;
k_121 = -(1/2)*(-2*zeta*rho*n+2*zeta*rho*n*alpha-1-2*rho*n*alpha)/(1+zeta*rho*n+rho*n);
k = k_121;
Efficiency_0_matrix = [k^2*xi/(2*zeta*rho*n*(k - 1 + alpha) + xi*k^2), 2*n*rho*(k - 1 + alpha)*zeta/(2*zeta*rho*n*(k - 1 + alpha) + xi*k^2), (k - 1)^2/(k^2 + (-2*n*rho - 2)*k + 2*rho*n*alpha + 1), -2*(-alpha + k)*n*rho/(k^2 + (-2*n*rho - 2)*k + 2*rho*n*alpha + 1)];

Stage_0(1,7) =Efficiency_0_matrix(1);  %Concrete Compresion
Stage_0(1,8) =Efficiency_0_matrix(2);  %Rebar Compresion 
Stage_0(1,9) =Efficiency_0_matrix(3);  %Concrete Tension
Stage_0(1,10) =Efficiency_0_matrix(4); %Rebar Tension
Stage_0(1,11) =sum(Efficiency_0_matrix); %for receck


%........................
% shared stage 121
%........................
beta_121 = 1;
k_121 = -(1/2)*(-2*zeta*rho*n+2*zeta*rho*n*alpha-1-2*rho*n*alpha)/(1+zeta*rho*n+rho*n);
k = k_121;
m_121 = -(2*((3+3*rho*n+3*zeta*rho*n)*k^2+(-6*rho*n*alpha+6*zeta*rho*n*alpha-3-6*zeta*rho*n)*k+3*zeta*rho*n+1+3*rho*n*alpha^2-6*zeta*rho*n*alpha+3*zeta*rho*n*alpha^2))/(k-1);
%Efficincy_1_matrix := [eff_con_c_1, eff_rebar_c_1, eff_con_t_1, eff_rebar_t_1]
Efficiency_121_matrix = [k^2*xi/(2*zeta*rho*n*(k - 1 + alpha) + xi*k^2), 2*n*rho*(k - 1 + alpha)*zeta/(2*zeta*rho*n*(k - 1 + alpha) + xi*k^2), (k - 1)^2/(k^2 + (-2*n*rho - 2)*k + 2*rho*n*alpha + 1), -2*(-alpha + k)*n*rho/(k^2 + (-2*n*rho - 2)*k + 2*rho*n*alpha + 1)];
Stage_121(1,1)= beta_121;
Stage_121(1,2)= k_121 ;
Stage_121(1,3) = m_121 ;
Stage_121(1,4) =(1/2)*beta_121/(1-k_121);
Stage_121(1,5) = m_121*Mcr ;
Stage_121(1,6) =((1/2)*beta_121/(1-k_121))*phicr;

Stage_121(1,7) =Efficiency_121_matrix(1);  %Concrete Compresion
Stage_121(1,8) =Efficiency_121_matrix(2);  %Rebar Compresion 
Stage_121(1,9) =Efficiency_121_matrix(3);  %Concrete Tension
Stage_121(1,10) =Efficiency_121_matrix(4); %Rebar Tension
Stage_121(1,11) =sum(Efficiency_121_matrix); %for receck

%........................
% shared stage 2122
%........................
beta_2122 =(sqrt((alpha - 1)^2*(n^2*kappa^2*(zeta + 1)^2*rho^2 + 2*n*(((zeta - 1)*(eta - 1)*alpha + eta + zeta)*kappa^2 - 3*(eta - 1)*((zeta - 1/3)*alpha - zeta/3 + 1/3)*kappa + 2*(eta - 1)*zeta*(alpha - 1/2))*rho + eta*kappa^2 - 2*alpha*(eta - 1)*kappa + (-2 + 2*eta)*alpha - eta + 1)) + (-eta + 1)*alpha^2 + (2*eta + (3*zeta - 1)*kappa*n*rho + kappa - 2)*alpha - eta + 1 + (-zeta + 1)*kappa*n*rho)/((-eta + 1)*alpha^2 + (4*n*rho*zeta + 2*eta)*alpha - 2*n*rho*zeta - eta) ;
beta = beta_2122;
k_2122  = (alpha*beta - kappa)/(beta - kappa);
k = k_2122 ;
m_2122 =(((-2 + 2*eta)*k^3 + (-6*n*rho*zeta - 6*eta)*k^2 + (6*eta + (-12*alpha + 12)*zeta*n*rho)*k - 2*eta + (-6*alpha^2 + 12*alpha - 6)*zeta*n*rho)*beta^3 - 3*((eta - 1)*k^2 + (2*kappa*n*rho - 2*eta + 2)*k - 2*alpha*kappa*n*rho + eta - 1)*(k - 1)*beta^2 + (k - 1)^3*(eta - 1))/((k - 1)*beta^2) ;
Efficiency_2122_matrix = [k^2*xi/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), 2*n*rho*(k - 1 + alpha)*zeta/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), (k - 1)^2*(beta^2*eta + (-2*eta + 2)*beta + eta - 1)/((eta*(k - 1)^2 - 2*n*rho*(-alpha + k))*beta^2 - 2*(k - 1)^2*(eta - 1)*beta + (k - 1)^2*(eta - 1)), -2*rho*n*(-alpha + k)*beta^2/((eta*k^2 + (-2*n*rho - 2*eta)*k + 2*n*rho*alpha + eta)*beta^2 - 2*(k - 1)^2*(eta - 1)*beta + (k - 1)^2*(eta - 1))];
Stage_2122(1,1)= beta_2122 ;
Stage_2122(1,2)= k_2122  ;
Stage_2122(1,3) = m_2122 ;
Stage_2122(1,4) =(1/2)*beta_2122/(1-k_2122);
Stage_2122(1,5) = m_2122*Mcr ;
Stage_2122(1,6) =((1/2)*beta_2122/(1-k_2122))*phicr;

Stage_2122(1,7) =Efficiency_2122_matrix(1);  %Concrete Compresion
Stage_2122(1,8) =Efficiency_2122_matrix(2);  %Rebar Compresion 
Stage_2122(1,9) =Efficiency_2122_matrix(3);  %Concrete Tension
Stage_2122(1,10) =Efficiency_2122_matrix(4); %Rebar Tension
Stage_2122(1,11) =sum(Efficiency_2122_matrix); %for receck


%........................
% shared stage 2242
%........................
beta_2242 = tau;
beta = beta_2242;
k_2242  = (-sqrt((2*(n*rho*zeta*alpha + 1/2)*(-1 + tau)^2*eta + n^2*(tau*zeta + kappa)^2*rho^2 - 2*(((-1 + tau)^2*alpha - tau^2)*zeta - tau*kappa)*n*rho + 2*tau - 1)*tau^2) + (n*rho*zeta + eta)*tau^2 + (kappa*n*rho - 2*eta + 2)*tau + eta - 1)/((-1 + tau)^2*(eta - 1));
k = k_2242 ;
m_2242 = (((-2 + 2*eta)*k^3 + (-6*n*rho*zeta - 6*eta)*k^2 + (6*eta + (-12*alpha + 12)*zeta*n*rho)*k - 2*eta + (-6*alpha^2 + 12*alpha - 6)*zeta*n*rho)*beta^3 - 3*((eta - 1)*k^2 + (2*kappa*n*rho - 2*eta + 2)*k - 2*alpha*kappa*n*rho + eta - 1)*(k - 1)*beta^2 + (k - 1)^3*(eta - 1))/((k - 1)*beta^2);
Efficiency_2242_matrix =[k^2*xi/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), 2*n*rho*(k - 1 + alpha)*zeta/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), (beta^2*eta + (-2*eta + 2)*beta + eta - 1)*(k - 1)/(eta*(k - 1)*beta^2 + ((-2*k + 2)*eta - 2*n*rho*kappa + 2*k - 2)*beta + (eta - 1)*(k - 1)), -2*rho*n*kappa*beta/(eta*(k - 1)*beta^2 + ((-2*k + 2)*eta - 2*n*rho*kappa + 2*k - 2)*beta + (eta - 1)*(k - 1))];
Stage_2242(1,1)= beta_2242 ;
Stage_2242(1,2)= k_2242  ;
Stage_2242(1,3) = m_2242 ;
Stage_2242(1,4) =(1/2)*beta_2242/(1-k_2242);
Stage_2242(1,5) = m_2242*Mcr ;
Stage_2242(1,6) =((1/2)*beta_2242/(1-k_2242))*phicr;

Stage_2242(1,7) =Efficiency_2242_matrix(1);  %Concrete Compresion
Stage_2242(1,8) =Efficiency_2242_matrix(2);  %Rebar Compresion 
Stage_2242(1,9) =Efficiency_2242_matrix(3);  %Concrete Tension
Stage_2242(1,10) =Efficiency_2242_matrix(4); %Rebar Tension
Stage_2242(1,11) =sum(Efficiency_2242_matrix); %for receck

%........................
% shared stage 2141
%........................
beta_2141 = tau;
k_2141 = -(-zeta*rho*n*tau^2+2*eta*tau-2*tau-eta+1-rho*n*tau^2-tau^2*eta+sqrt(rho^2*n^2*tau^4-2*rho*n*tau^2+tau^2*eta+tau^4*eta-2*tau^3*eta+4*rho*n*tau^3-2*eta*rho*n*tau^2*alpha+4*zeta*rho*n*tau^3*alpha+4*rho*n*tau^3*alpha*eta-2*rho*n*tau^4*alpha*eta-2*zeta*rho*n*tau^2*alpha-2*zeta*rho*n*tau^4*alpha+2*zeta*rho*n*tau^4*alpha*eta-4*zeta*rho*n*tau^3*alpha*eta+2*eta*zeta*rho*n*tau^2*alpha+2*tau^3-tau^2-4*rho*n*tau^3*alpha+2*rho*n*tau^4*alpha+2*tau^4*zeta*rho*n+2*eta*rho*n*tau^2+2*rho*n*tau^4*eta-4*rho*n*tau^3*eta+2*rho*n*tau^2*alpha+zeta^2*rho^2*n^2*tau^4+2*zeta*rho^2*n^2*tau^4))/(-2*eta*tau-tau^2+eta-1+2*tau+tau^2*eta);
k=k_2141;
mm_2141 = (((-2+2*eta)*tau^3-(3*(eta-1))*tau^2+eta-1)*k^3+((-6*eta+(-6*zeta-6)*n*rho)*tau^3-(3*(-3*eta+3))*tau^2-3*eta+3)*k^2+((6*eta+((-12*alpha+12)*zeta+12*alpha)*n*rho)*tau^3-(3*(3*eta-3))*tau^2+3*eta-3)*k+(-2*eta+((-6*alpha^2+12*alpha-6)*zeta-6*alpha^2)*n*rho)*tau^3-(3*(-eta+1))*tau^2-eta+1)/((k-1)*tau^2);
Efficiency_2141_matrix = [k^2*xi/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), 2*n*rho*(k - 1 + alpha)*zeta/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), (k - 1)^2*(beta^2*eta + (-2*eta + 2)*beta + eta - 1)/((eta*(k - 1)^2 - 2*n*rho*(-alpha + k))*beta^2 - 2*(k - 1)^2*(eta - 1)*beta + (k - 1)^2*(eta - 1)), -2*rho*n*(-alpha + k)*beta^2/((eta*k^2 + (-2*n*rho - 2*eta)*k + 2*n*rho*alpha + eta)*beta^2 - 2*(k - 1)^2*(eta - 1)*beta + (k - 1)^2*(eta - 1))];
Stage_2141(1,1) = beta_2141;
Stage_2141(1,2) = k_2141 ;
Stage_2141(1,3) = mm_2141;
Stage_2141(1,4) =(1/2)*beta_2141/(1-k_2141);
Stage_2141(1,5) = mm_2141*Mcr ;
Stage_2141(1,6) =((1/2)*beta_2141/(1-k_2141))*phicr;

Stage_2141(1,7) =Efficiency_2141_matrix(1);  %Concrete Compresion
Stage_2141(1,8) =Efficiency_2141_matrix(2);  %Rebar Compresion 
Stage_2141(1,9) =Efficiency_2141_matrix(3);  %Concrete Tension
Stage_2141(1,10) =Efficiency_2141_matrix(4); %Rebar Tension
Stage_2141(1,11) =sum(Efficiency_2141_matrix); %for receck
%........................

%........................
% Set1 4142, 4252
%........................
beta_4142=(mu*alpha^2+rho*n*kappa-2*alpha*mu+3*kappa*zeta*rho*n*alpha+mu+kappa*alpha-kappa*zeta*rho*n-kappa*rho*n*alpha+sqrt((alpha-1)^2*(6*rho*n*alpha*mu*kappa*zeta-8*eta*tau*zeta*rho*n*alpha-8*mu*tau*zeta*rho*n*alpha+4*eta*tau^2*zeta*rho*n*alpha+n^2*rho^2*kappa^2*zeta^2+2*n^2*rho^2*zeta*kappa^2-2*eta*zeta*rho*n-4*tau*zeta*rho*n-4*zeta*rho*n*alpha+2*rho*n*alpha*kappa^2+2*zeta*rho*n*kappa^2-2*mu^2*alpha-alpha^2+mu^2-2*rho*n*alpha*mu*kappa+4*eta*tau*zeta*rho*n+4*mu*tau*zeta*rho*n-2*eta*tau^2*zeta*rho*n-2*rho*n*kappa*mu*zeta+8*tau*zeta*rho*n*alpha+4*eta*zeta*rho*n*alpha-2*zeta*rho*n*alpha*kappa^2+alpha^2*mu^2+eta*alpha^2+2*tau*alpha^2+2*rho*n*kappa*mu-2*mu*tau*alpha^2-2*eta*tau*alpha^2+eta*tau^2*alpha^2+2*zeta*rho*n+2*mu*kappa*alpha+n^2*rho^2*kappa^2)))/(-2*zeta*rho*n+alpha^2+4*zeta*rho*n*alpha);
beta = beta_4142;
k_4142 =(beta*alpha-kappa)/(beta-kappa);
k=k_4142;
mm_4142 =((2*eta*tau^3+(-3*mu-3*eta+3)*tau^2-2*beta^3*xi+3*beta^2*mu+eta-1)*k^3+(-6*eta*tau^3+(9*mu+9*eta-9)*tau^2-6*n*rho*(zeta+1)*beta^3-9*beta^2*mu-3*eta+3)*k^2+(6*eta*tau^3+(-9*mu-9*eta+9)*tau^2-12*n*((alpha-1)*zeta-alpha)*rho*beta^3+9*beta^2*mu+3*eta-3)*k-2*eta*tau^3+(3*mu+3*eta-3)*tau^2-(6*((alpha-1)^2*zeta+alpha^2))*n*rho*beta^3-3*beta^2*mu-eta+1)/(beta^2*(k-1));
Efficiency_4142_matrix = [k^2*xi/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), 2*n*rho*(k - 1 + alpha)*zeta/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), 2*(k - 1)^2*(-eta*tau^2/2 + (mu + eta - 1)*tau - mu*beta - eta/2 + 1/2)/(-eta*(k - 1)^2*tau^2 + 2*(k - 1)^2*(mu + eta - 1)*tau + (-2*beta*mu - eta + 1)*k^2 + (2*beta^2*n*rho + 4*beta*mu + 2*eta - 2)*k - 2*rho*n*beta^2*alpha - 2*mu*beta - eta + 1), 2*rho*n*(-alpha + k)*beta^2/((-eta*tau^2 + (2*mu + 2*eta - 2)*tau - 2*mu*beta - eta + 1)*k^2 + (2*eta*tau^2 + (-4*mu - 4*eta + 4)*tau + 2*n*rho*beta^2 + 4*mu*beta + 2*eta - 2)*k - eta*tau^2 + (2*mu + 2*eta - 2)*tau - 2*rho*n*beta^2*alpha - 2*mu*beta - eta + 1)];
Stage_4142(1,1) = beta_4142;
Stage_4142(1,2) = k_4142 ;
Stage_4142(1,3) = mm_4142;
Stage_4142(1,4) =(1/2)*beta_4142/(1-k_4142);
Stage_4142(1,5) = mm_4142*Mcr ;
Stage_4142(1,6) =((1/2)*beta_4142/(1-k_4142))*phicr;

Stage_4142(1,7) =Efficiency_4142_matrix(1);  %Concrete Compresion
Stage_4142(1,8) =Efficiency_4142_matrix(2);  %Rebar Compresion 
Stage_4142(1,9) =Efficiency_4142_matrix(3);  %Concrete Tension
Stage_4142(1,10) =Efficiency_4142_matrix(4); %Rebar Tension
Stage_4142(1,11) =sum(Efficiency_4142_matrix); %for receck
%........................

%........................
beta_4252  = -(1/2)*(-omega*zeta*rho*n-mu+2*zeta*rho*n*alpha*omega-rho*n*kappa+sqrt(2*zeta*rho*n+rho^2*n^2*kappa^2+rho^2*n^2*zeta^2*omega^2-4*eta*tau*zeta*rho*n*alpha-4*mu*tau*zeta*rho*n*alpha+2*eta*tau^2*zeta*rho*n*alpha-2*zeta*rho*n*alpha+2*omega*zeta*rho*n*mu-2*omega^2*zeta*rho*n*alpha+4*tau*zeta*rho*n*alpha+2*eta*zeta*rho*n*alpha-2*rho^2*n^2*zeta*kappa*omega-4*tau*zeta*rho*n-4*zeta*rho*n*alpha*omega*mu+2*omega^2*zeta*rho*n+4*eta*tau*zeta*rho*n-2*eta*tau^2*zeta*rho*n+2*rho*n*kappa*mu+mu^2+4*mu*tau*zeta*rho*n-2*eta*zeta*rho*n))/(zeta*rho*n*(alpha-1));
beta = beta_4252;
k_4252 =omega/(omega+beta);
k=k_4252;
mm_4252 = ((3*tau^2+3*mu*beta^2+2*eta*tau^3+eta-3*mu*tau^2-3*eta*tau^2-1-2*xi*beta^3)*k^3+(3+9*mu*tau^2+9*eta*tau^2-6*eta*tau^3-6*zeta*rho*n*beta^3-3*eta-6*rho*n*kappa*beta^2-9*tau^2-9*mu*beta^2)*k^2+(-12*zeta*rho*n*beta^3*alpha+9*mu*beta^2+9*tau^2-9*mu*tau^2+6*rho*n*kappa*beta^2*alpha+6*rho*n*kappa*beta^2+3*eta-9*eta*tau^2+12*zeta*rho*n*beta^3+6*eta*tau^3-3)*k+1-eta-3*tau^2+3*mu*tau^2-6*zeta*rho*n*beta^3-6*rho*n*kappa*beta^2*alpha-3*mu*beta^2-6*zeta*rho*n*beta^3*alpha^2-2*eta*tau^3+3*eta*tau^2+12*zeta*rho*n*beta^3*alpha)/(beta^2*(k-1));
Efficiency_4252_matrix =[k^2*xi/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), 2*n*rho*(k - 1 + alpha)*zeta/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), 2*(-eta*tau^2/2 + (mu + eta - 1)*tau - mu*beta - eta/2 + 1/2)*(k - 1)/(-eta*(k - 1)*tau^2 + 2*(mu + eta - 1)*(k - 1)*tau + (1 - k)*eta + (-2*beta*mu + 1)*k + 2*beta*kappa*n*rho + 2*mu*beta - 1), 2*n*rho*beta*kappa/(-eta*(k - 1)*tau^2 + 2*(mu + eta - 1)*(k - 1)*tau + (-2*beta*mu - eta + 1)*k + 2*beta*kappa*n*rho + 2*mu*beta + eta - 1)];
Stage_4252(1,1) = beta_4252;
Stage_4252(1,2) = k_4252;
Stage_4252(1,3) = mm_4252;
Stage_4252(1,4) =(1/2)*beta_4252/(1-k_4252);
Stage_4252(1,5) = mm_4252*Mcr ;
Stage_4252(1,6) =((1/2)*beta_4252/(1-k_4252))*phicr;

Stage_4252(1,7) =Efficiency_4252_matrix(1);  %Concrete Compresion
Stage_4252(1,8) =Efficiency_4252_matrix(2);  %Rebar Compresion 
Stage_4252(1,9) =Efficiency_4252_matrix(3);  %Concrete Tension
Stage_4252(1,10) =Efficiency_4252_matrix(4); %Rebar Tension
Stage_4252(1,11) =sum(Efficiency_4252_matrix); %for receck
%.......end.............

%........................
% Set1 4151, 5152
%........................
beta_4151=-(1/2)*(-mu+rho*n*omega+2*zeta*rho*n*alpha*omega-2*rho*n*alpha*omega-zeta*rho*n*omega+sqrt(-4*zeta*rho*n*alpha*omega*mu-4*eta*n*rho*zeta*alpha*tau+2*tau^2*zeta*rho*n*alpha*eta-4*mu*tau*zeta*rho*n*alpha+4*mu*n*rho*alpha*omega+4*eta*n*rho*tau*zeta-2*eta*n*rho*tau^2*zeta+4*mu*tau*zeta*rho*n+2*mu*n*rho*zeta*omega-2*eta*n*rho*tau^2*alpha+4*eta*n*rho*tau*alpha+4*rho*n*tau*zeta*alpha+4*mu*tau*rho*n*alpha+2*eta*n*rho*zeta*alpha-2*zeta*rho*n*omega^2*alpha+2*rho*n*alpha+omega^2*rho^2*n^2+2*zeta*rho*n+mu^2-4*rho*n*tau*zeta-2*eta*n*rho*zeta+2*zeta*rho*n*omega^2+2*n^2*rho^2*zeta*omega^2-2*mu*n*rho*omega-4*rho*n*alpha*tau-2*rho*n*alpha*eta+2*omega^2*rho*n*alpha-2*zeta*rho*n*alpha+omega^2*rho^2*n^2*zeta^2))/(rho*n*(-alpha+zeta*alpha-zeta));
beta = beta_4151;
k_4151 =omega/(omega+beta);
k = k_4151;
mm_4151 = ((2*eta*tau^3+(-3*mu-3*eta+3)*tau^2-2*beta^3*xi+3*beta^2*mu+eta-1)*k^3+(-6*eta*tau^3+(9*mu+9*eta-9)*tau^2-6*n*rho*(zeta+1)*beta^3-9*beta^2*mu-3*eta+3)*k^2+(6*eta*tau^3+(-9*mu-9*eta+9)*tau^2-12*n*((alpha-1)*zeta-alpha)*rho*beta^3+9*beta^2*mu+3*eta-3)*k-2*eta*tau^3+(3*mu+3*eta-3)*tau^2-(6*((alpha-1)^2*zeta+alpha^2))*n*rho*beta^3-3*beta^2*mu-eta+1)/(beta^2*(k-1));
Efficiency_4151_matrix = [k^2*xi/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), 2*n*rho*(k - 1 + alpha)*zeta/(2*zeta*rho*n*(k - 1 + alpha) + k^2*xi), 2*(k - 1)^2*(-eta*tau^2/2 + (mu + eta - 1)*tau - mu*beta - eta/2 + 1/2)/(-eta*(k - 1)^2*tau^2 + 2*(k - 1)^2*(mu + eta - 1)*tau + (-2*beta*mu - eta + 1)*k^2 + (2*beta^2*n*rho + 4*beta*mu + 2*eta - 2)*k - 2*rho*n*beta^2*alpha - 2*mu*beta - eta + 1), 2*rho*n*(-alpha + k)*beta^2/((-eta*tau^2 + (2*mu + 2*eta - 2)*tau - 2*mu*beta - eta + 1)*k^2 + (2*eta*tau^2 + (-4*mu - 4*eta + 4)*tau + 2*n*rho*beta^2 + 4*mu*beta + 2*eta - 2)*k - eta*tau^2 + (2*mu + 2*eta - 2)*tau - 2*rho*n*beta^2*alpha - 2*mu*beta - eta + 1)];
Stage_4151(1,1) = beta_4151;
Stage_4151(1,2) = k_4151;
Stage_4151(1,3) = mm_4151;
Stage_4151(1,4) =(1/2)*beta_4151/(1-k_4151);
Stage_4151(1,5) = mm_4151*Mcr ;
Stage_4151(1,6) =((1/2)*beta_4151/(1-k_4151))*phicr;

Stage_4151(1,7) =Efficiency_4151_matrix(1);  %Concrete Compresion
Stage_4151(1,8) =Efficiency_4151_matrix(2);  %Rebar Compresion 
Stage_4151(1,9) =Efficiency_4151_matrix(3);  %Concrete Tension
Stage_4151(1,10) =Efficiency_4151_matrix(4); %Rebar Tension
Stage_4151(1,11) =sum(Efficiency_4151_matrix); %for receck
%.......end.............




beta_5152 = (1/2)*(n*rho*kappa-2*mu*alpha+3*alpha*zeta*n*rho*kappa-alpha*n*rho*kappa+mu-alpha*omega-zeta*rho*n*kappa+alpha^2*omega+mu*alpha^2+sqrt((-1+alpha)^2*(4*alpha*zeta*n*rho*eta*tau^2-8*alpha*zeta*n*rho*mu*tau-8*alpha*zeta*n*rho*eta*tau+6*alpha*zeta*n*rho*kappa*omega+6*alpha*zeta*n*rho*kappa*mu+2*alpha^2*omega*mu-2*alpha*omega*mu+2*zeta*rho*n+rho^2*n^2*kappa^2-4*alpha*zeta*n*rho+2*n*rho*kappa*mu-4*zeta*rho*n*tau-2*zeta*rho*n*omega^2-2*zeta*rho*n*eta+rho^2*n^2*zeta^2*kappa^2+2*zeta*rho^2*n^2*kappa^2+mu^2+4*alpha*zeta*n*rho*eta-2*zeta*rho*n*eta*tau^2+4*zeta*rho*n*mu*tau+4*zeta*rho*n*eta*tau-4*zeta*rho*n*kappa*omega-2*zeta*rho*n*kappa*mu-2*alpha*n*rho*kappa*omega-2*alpha*n*rho*kappa*mu+8*alpha*zeta*n*rho*tau+4*alpha*zeta*n*rho*omega^2+alpha^2*mu^2+omega^2*alpha^2-2*alpha*mu^2)))/(zeta*rho*n*(2*alpha-1));
beta =beta_5152;
k_5152=(beta*alpha-kappa)/(beta-kappa);
k=k_5152;
mm_5152=((2*eta*tau^3-3*eta*tau^2+eta+3*tau^2-1+3*mu*beta^2-3*mu*tau^2+3*omega*beta^2-omega^3)*k^3+(3-9*tau^2+9*mu*tau^2-3*omega*beta^2-6*zeta*rho*n*beta^3+9*eta*tau^2-9*mu*beta^2+3*omega^3-6*rho*n*kappa*beta^2-3*eta-6*eta*tau^3)*k^2+(9*tau^2-9*eta*tau^2+6*rho*n*kappa*beta^2*alpha-3-12*zeta*rho*n*beta^3*alpha+6*eta*tau^3+9*mu*beta^2-3*omega^3+6*rho*n*kappa*beta^2-9*mu*tau^2+3*eta+12*zeta*rho*n*beta^3)*k+1-eta-6*rho*n*kappa*beta^2*alpha-3*mu*beta^2+3*mu*tau^2-2*eta*tau^3+3*eta*tau^2+omega^3-3*tau^2-6*zeta*rho*n*beta^3*alpha^2-6*zeta*rho*n*beta^3+12*zeta*rho*n*beta^3*alpha)/(beta^2*(-1+k));
Efficiency_5152_matrix = [(xi*omega^2*(k - 1)^2 + 2*omega*k*xi*(k - 1)*beta)/(xi*omega^2*(k - 1)^2 + 2*omega*k*xi*(k - 1)*beta - 2*zeta*rho*n*(k - 1 + alpha)*beta^2), -2*n*rho*beta^2*zeta*(k - 1 + alpha)/(xi*omega^2*(k - 1)^2 + 2*omega*k*xi*(k - 1)*beta - 2*zeta*rho*n*(k - 1 + alpha)*beta^2), 2*(k - 1)^2*(-eta*tau^2/2 + (mu + eta - 1)*tau - mu*beta - eta/2 + 1/2)/(-eta*(k - 1)^2*tau^2 + 2*(k - 1)^2*(mu + eta - 1)*tau + (-2*beta*mu - eta + 1)*k^2 + (2*beta^2*n*rho + 4*beta*mu + 2*eta - 2)*k - 2*rho*n*beta^2*alpha - 2*mu*beta - eta + 1), 2*rho*n*(-alpha + k)*beta^2/((-eta*tau^2 + (2*mu + 2*eta - 2)*tau - 2*mu*beta - eta + 1)*k^2 + (2*eta*tau^2 + (-4*mu - 4*eta + 4)*tau + 2*n*rho*beta^2 + 4*mu*beta + 2*eta - 2)*k - eta*tau^2 + (2*mu + 2*eta - 2)*tau - 2*rho*n*beta^2*alpha - 2*mu*beta - eta + 1)];
Stage_5152(1,1) = beta_5152;
Stage_5152(1,2) = k_5152;
Stage_5152(1,3) = mm_5152;
Stage_5152(1,4) =(1/2)*beta_5152/(1-k_5152);
Stage_5152(1,5) = mm_5152*Mcr ;
Stage_5152(1,6) =((1/2)*beta_5152/(1-k_5152))*phicr;

Stage_5152(1,7) =Efficiency_5152_matrix(1);  %Concrete Compresion
Stage_5152(1,8) =Efficiency_5152_matrix(2);  %Rebar Compresion 
Stage_5152(1,9) =Efficiency_5152_matrix(3);  %Concrete Tension
Stage_5152(1,10) =Efficiency_5152_matrix(4); %Rebar Tension
Stage_5152(1,11) =sum(Efficiency_5152_matrix); %for receck
%.......end.............


%........................
% Infinity
%........................
k_inf = (-zeta*rho*n*kappa+rho*n*kappa+mu)/(omega+mu);
k=k_inf;
beta=100000000000000000000000000000000000000;
m_inf = (3*(-2*rho*n*kappa*mu+omega*mu+2*rho*n*kappa*mu*alpha+2*mu*zeta*rho*n*kappa*alpha-2*omega*zeta*rho*n*kappa+2*rho*n*kappa*alpha*omega+2*zeta*rho*n*kappa*alpha*omega+2*zeta*rho^2*n^2*kappa^2-rho^2*n^2*kappa^2-zeta^2*rho^2*n^2*kappa^2))/(omega+mu);
Efficiency_53inf_matrix = [(xi*omega^2*(k - 1) + 2*omega*k*xi*beta)/(xi*omega^2*(k - 1) + 2*omega*k*xi*beta + 2*zeta*rho*n*kappa*beta), 2*n*rho*zeta*beta*kappa/(xi*omega^2*(k - 1) + 2*omega*k*xi*beta + 2*zeta*rho*n*kappa*beta), 2*(-eta*tau^2/2 + (mu + eta - 1)*tau - mu*beta - eta/2 + 1/2)*(k - 1)/(-eta*(k - 1)*tau^2 + 2*(mu + eta - 1)*(k - 1)*tau + (1 - k)*eta + (-2*beta*mu + 1)*k + 2*beta*kappa*n*rho + 2*mu*beta - 1), 2*n*rho*beta*kappa/(-eta*(k - 1)*tau^2 + 2*(mu + eta - 1)*(k - 1)*tau + (-2*beta*mu - eta + 1)*k + 2*beta*kappa*n*rho + 2*mu*beta + eta - 1)];
Stage_inf(1,1) = 0;
Stage_inf(1,2) = k_inf ;
Stage_inf(1,3) = m_inf;
Stage_inf(1,4) =0;
Stage_inf(1,5) = m_inf*Mcr ;
Stage_inf(1,6) =0*phicr;

Stage_inf(1,7) =Efficiency_53inf_matrix(1);  %Concrete Compresion
Stage_inf(1,8) =Efficiency_53inf_matrix(2);  %Rebar Compresion 
Stage_inf(1,9) =Efficiency_53inf_matrix(3);  %Concrete Tension
Stage_inf(1,10) =Efficiency_53inf_matrix(4); %Rebar Tension
Stage_inf(1,11) =sum(Efficiency_53inf_matrix); %for receck
%.......end.............




%-------------------------------------------------------------------------
%  Output File
%-------------------------------------------------------------------------

Output(1:10,1) = [0,121,2122,2242,2141,4142,4252,4151,5152,inf].' ;
Output(1,2:12) = Stage_0;
Output(2,2:12) = Stage_121;
Output(3,2:12) = Stage_2122;
Output(4,2:12) = Stage_2242;
Output(5,2:12) = Stage_2141;
Output(6,2:12) = Stage_4142;
Output(7,2:12) = Stage_4252;
Output(8,2:12) = Stage_4151;
Output(9,2:12) = Stage_5152;
Output(10,2:12) = Stage_inf;


Output_dummy_case1 = [Stage_0,Stage_121,Stage_2141,Stage_4151,Stage_5152];
Output_dummy_case2 = [Stage_0,Stage_121,Stage_4142,Stage_4251,Stage_4252];
Output_dummy_case3 = [Stage_0,Stage_121,Stage_2122,Stage_2242,Stage_4252];



if Output(6,2) >= Output(8,2)% || Output(6,2) <= Output(5,2)
   Output_Final(1:5,1)  = [0,121,2141,4151,5152].';
   for i=1:5
       Output_Final(i,2:12) = Output_dummy_case1(i*11-10:i*11);
       casee = 1;
   end
elseif  Output(3,2) <= Output(5,2) 
   Output_Final(1:5,1)  = [0,121,2122,2242,4252].';
   for k=1:5
       Output_Final(k,2:12) = Output_dummy_case3(k*11-10:k*11);
       casee = 3;
   end
else 
   Output_Final(1:5,1)  = [0,121,2141,4142,4252].';
   for j=1:5
       Output_Final(j,2:12) = Output_dummy_case2(j*11-10:j*11);
       casee = 2;
   end
end 

% Replace data that is NAN with Zero
Output_Final(isnan( Output_Final)) = 9999999999;