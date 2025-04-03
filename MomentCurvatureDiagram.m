function btldkpnfcvmm = MomentCurvatureDiagram(ecr,E,alpha,beta_tu,mu,eta,gamma,omega,lambda_cu,subMC)
% Got fomular for kp,nf,cv,mm from Bilinear_Residual_FRC_model_02.mw

ecu = lambda_cu*ecr;
etu = beta_tu*ecr;

b1  = linspace(0,1,subMC(1));
b2  = linspace(1,alpha,subMC(2)+1);
b3  = linspace(alpha,beta_tu,subMC(3)+1);
bt  = [b1(1:end),b2(2:end),b3(2:end)];       % remove a duplicate node at joining region 1-2 and 2-3

nInc= length(bt);
ld  = zeros(1,nInc);
kp  = zeros(1,nInc);
nf  = zeros(1,nInc);
mm  = zeros(1,nInc);
cv  = zeros(1,nInc);

for i=2:nInc
    beta = bt(i);
    if beta<=1
        if gamma == 1
            kp(i) = 0.5;
        else
            kp(i) = (-1+gamma^(1/2))/(-1+gamma);
        end
        k     = kp(i);
        nf(i) = 1/2*beta*(-1+(-1+gamma)*k^2+2*k)/(-1+k);
        mm(i) = -2*beta*(1+(-1+gamma)*k^3+3*k^2-3*k)/(-1+k);
        cv(i) = -1/2*beta/(-1+k);
        
    elseif beta<=alpha
        K21   = -1+2*beta+eta*beta^2-2*eta*beta+eta-beta^2*gamma;
        kp(i) = (beta^2*gamma+K21-(gamma^2*beta^4+K21*gamma*beta^2)^(1/2))/K21;
        k     = kp(i);
        lambda = beta*k/(1-k);
        
        if lambda <= omega
            NN21  =  (-1+eta*beta^2+2*beta+eta-2*eta*beta)/beta;
            MM21  =  (-3*beta^2-eta-2*eta*beta^3+3*eta*beta^2+1)/beta^2;
            nf(i) = -(k^2*(NN21-beta*gamma)-2*k*NN21+NN21)/(-2+2*k);
            mm(i) = -((2*beta*gamma+MM21)*k^3-3*k^2*MM21+3*k*MM21-MM21)/(-1+k);
            cv(i) = -1/2*beta/(-1+k);   
        else
            K22   = -1+2*beta+eta*beta^2-2*eta*beta+omega^2*gamma+eta;
            kp(i) = K22/(K22+2*omega*gamma*beta);
            k     = kp(i);
            lambda = beta*k/(1-k);
            if lambda <= lambda_cu
                NN22  = (1-2*beta-eta*beta^2+2*eta*beta-omega^2*gamma-eta)/(2*beta);
                MM22  = (eta-1+3*beta^2-omega^3*gamma-3*eta*beta^2+2*eta*beta^3)/beta^2;
                nf(i) = (NN22-omega*gamma)*k-NN22;
                mm(i) = (3*omega*gamma+MM22)*k^2-2*k*MM22+MM22;
                cv(i) = -1/2*beta/(-1+k);
            else
                nf(i) = 0.0;
                mm(i) = 0.0;
                cv(i) = (1+10^-3)*cv(i-1);
            end
        end
        
    elseif beta<=beta_tu
        K31   = -2*mu*alpha-2*eta*alpha+eta+eta*alpha^2+2*alpha-1+2*mu*beta;
        kp(i) = (K31-(gamma*K31*beta^2)^(1/2))/(-beta^2*gamma+K31);
        k     = kp(i);
        lambda = beta*k/(1-k);
        
        if lambda <= omega
            NN31  = (-eta*alpha^2-2*mu*beta+2*mu*alpha-eta+1-2*alpha+2*eta*alpha)/beta;
            MM31  = (-1+2*eta*alpha^3+3*alpha^2+3*mu*beta^2-3*eta*alpha^2+eta-3*mu*alpha^2)/beta^2;
            nf(i) = (k^2*(NN31+beta*gamma)-2*k*NN31+NN31)/(-2+2*k);
            mm(i) = ((MM31-2*beta*gamma)*k^3-3*MM31*k^2+3*k*MM31-MM31)/(-1+k);
            cv(i) = -1/2*beta/(-1+k);
        else
            K32   = -1+omega^2*gamma+2*alpha+eta*alpha^2-2*eta*alpha+eta+2*mu*beta-2*mu*alpha;
            kp(i) = K32/(K32+2*omega*gamma*beta);
            k     = kp(i);
            lambda = beta*k/(1-k);
            if lambda <= lambda_cu
                NN32  = (1-omega^2*gamma-2*alpha-eta*alpha^2+2*eta*alpha-eta-2*mu*beta+2*mu*alpha)/(2*beta);
                MM32  = (-1+eta-3*eta*alpha^2+2*eta*alpha^3+3*mu*beta^2-omega^3*gamma-3*mu*alpha^2+3*alpha^2)/beta^2;
                nf(i) =  k*(NN32-omega*gamma)-NN32;
                mm(i) = (MM32+3*omega*gamma)*k^2-2*k*MM32+MM32;
                cv(i) = -1/2*beta/(-1+k);
            else
                nf(i) = 0.0;
                mm(i) = 0.0;
                cv(i) = (1+10^-3)*cv(i-1);
            end
        end
    end
    ld(i) = beta*k/(1-k);    % take the final k to compute ld(i)
end
kp(1) = kp(2);

% trim data upto mm = 0
stopID = length(mm);
for i=2:length(mm)
    if mm(i) == 0
        stopID = i;
        break;
    end    
end

BT = bt(1:stopID);
LD = ld(1:stopID);
KP = kp(1:stopID);
NF = nf(1:stopID);
CV = cv(1:stopID);
MM = mm(1:stopID);

btldkpnfcvmm = [BT; LD; KP; NF; CV; MM]';




