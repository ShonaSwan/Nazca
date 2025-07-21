tic;

% define dimensional scales
h0    = D/10;
e0    = 1e18;
Drho0 = 500;
p0    = Drho0*g0*h0;
u0    = p0*h0./e0;
K0    = h0^2/e0;
t0    = h0/u0;


if ~bnchm && step>0 && ~restart

% residual of mixture mass evolution
% Div_rhoV = reshape(DDs*[W(:);U(:)] + DDm*[wm(:);um(:)] ,Nz+2,Nx+2) * Drho0/h0;
% drhodt   = - Div_rhoV(2:end-1,2:end-1) + Gem + Gex;

res_rho  = (a1*rho-a2*rhoo-a3*rhooo)/dt - (b1*drhodt + b2*drhodto + b3*drhodtoo);
% upd_rho  = - alpha*res_rho*dt/a1;

% upd_rho = - alpha*res_rho./b1;
% VolSrc  = VolSrc + upd_rho;
VolSrc  = - (a1*rho-a2*rhoo-a3*rhooo)/dt + Gem + Gex;  % correct volume source term by scaled residual

end


%% assemble coefficients for matrix velocity diagonal and right-hand side (KV and RV)

IIL = [];       % equation indeces into L
JJL = [];       % variable indeces into L
AAL = [];       % coefficients for L
IIR = [];       % equation indeces into R
AAR = [];       % forcing entries for R


% assemble coefficients of z-stress divergence

% top boundary
ii  = MapW(1,:); jj1 = ii;  
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
aa  = zeros(size(ii));
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii  = MapW(end,:); jj1 = ii; jj2 = MapW(end-1,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
aa  = zeros(size(ii)); 
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% left boundary
ii  = MapW((2:end-1),1); jj1 = ii; jj2 = MapW((2:end-1),2);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapW((2:end-1),end); jj1 = ii; jj2 = MapW((2:end-1),end-1);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];  AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];  AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];


% internal points
ii    = MapW(2:end-1,2:end-1);
EtaC1 = etaco(2:end-1,1:end-1)/e0;   EtaC2 = etaco(2:end-1,2:end)/e0;
EtaP1 = eta  (1:end-1,:      )/e0;   EtaP2 = eta  (2:end,:      )/e0;

% coefficients multiplying z-velocities W
%             top          ||         bottom          ||           left            ||          right
jj1 = MapW(1:end-2,2:end-1); jj2 = MapW(3:end,2:end-1); jj3 = MapW(2:end-1,1:end-2); jj4 = MapW(2:end-1,3:end);

aa  = 2/3*(EtaP1+EtaP2)/(h/h0)^2 + 1/2*(EtaC1+EtaC2)/(h/h0)^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % W on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-2/3*EtaP1(:)/(h/h0)^2];      % W one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-2/3*EtaP2(:)/(h/h0)^2];      % W one below
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-1/2*EtaC1(:)/(h/h0)^2];      % W one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-1/2*EtaC2(:)/(h/h0)^2];      % W one to the right

% coefficients multiplying x-velocities U
%         top left         ||        bottom left          ||       top right       ||       bottom right
jj1 = MapU(2:end-2,1:end-1); jj2 = MapU(3:end-1,1:end-1); jj3 = MapU(2:end-2,2:end); jj4 = MapU(3:end-1,2:end);

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP1(:))/(h/h0)^2];  % U one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;+(1/2*EtaC1(:)-1/3*EtaP2(:))/(h/h0)^2];  % U one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+(1/2*EtaC2(:)-1/3*EtaP1(:))/(h/h0)^2];  % U one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP2(:))/(h/h0)^2];  % U one to the bottom and right


% z-RHS vector
rr  = + (Drhow(2:end-1,:).*g0)/(Drho0*g0);
if bnchm; rr = rr + src_W_mms(2:end-1,:); end

IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];


%  assemble coefficients of x-stress divergence

% top boundary
ii  = MapU(1,:); jj1 = ii; jj2 = MapU(2,:);
aa  = zeros(size(ii)) + bnd_spr/u0 * 2;
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii  = MapU(end,:); jj1 = ii; jj2 = MapU(end-1,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% left boundary
ii  = MapU((2:end-1),1); jj = ii;
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa  = zeros(size(ii));% + UBG((2:end-1),1);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapU((2:end-1),end); jj1 = ii; jj2 = MapU((2:end-1),end-1);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
aa  = zeros(size(ii));
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% internal points
ii    = MapU(2:end-1,2:end-1);
EtaC1 = etaco(1:end-1,2:end-1)/e0;  EtaC2 = etaco(2:end,2:end-1)/e0;
EtaP1 = eta  (:      ,1:end-1)/e0;  EtaP2 = eta  (:      ,2:end)/e0;

% coefficients multiplying x-velocities U
%            left          ||          right          ||           top             ||          bottom
 jj1 = MapU(2:end-1,1:end-2); jj2 = MapU(2:end-1,3:end); jj3 = MapU(1:end-2,2:end-1); jj4 = MapU(3:end,2:end-1);

aa  = 2/3*(EtaP1+EtaP2)/(h/h0)^2 + 1/2*(EtaC1+EtaC2)/(h/h0)^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % U on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-2/3*EtaP1(:)/(h/h0)^2];      % U one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-2/3*EtaP2(:)/(h/h0)^2];      % U one to the right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-1/2*EtaC1(:)/(h/h0)^2];      % U one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-1/2*EtaC2(:)/(h/h0)^2];      % U one below

% coefficients multiplying z-velocities W
%         top left         ||        top right          ||       bottom left       ||       bottom right
 jj1 = MapW(1:end-1,2:end-2); jj2 = MapW(1:end-1,3:end-1); jj3 = MapW(2:end,2:end-2); jj4 = MapW(2:end,3:end-1);

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP1(:))/(h/h0)^2];  % W one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;+(1/2*EtaC1(:)-1/3*EtaP2(:))/(h/h0)^2];  % W one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+(1/2*EtaC2(:)-1/3*EtaP1(:))/(h/h0)^2];  % W one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP2(:))/(h/h0)^2];  % W one to the bottom and right


% x-RHS vector
rr = zeros(size(ii));
if bnchm
    rr = rr + src_U_mms(2:end-1,2:end-1);
end

IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];

% assemble coefficient matrix & right-hand side vector
KV  = sparse(IIL,JJL,AAL,NW+NU,NW+NU);
RV  = sparse(IIR,ones(size(IIR)),AAR);


%% assemble coefficients for gradient operator (GG)

if ~exist('GG','var') || bnchm
    IIL = [];       % equation indeces into A
    JJL = [];       % variable indeces into A
    AAL = [];       % coefficients for A
    
    
    % Coefficients for z-gradient
    ii  = MapW(2:end-1,2:end-1);
    
    %         top              ||          bottom
    jj1 = MapP(2:end-2,2:end-1); jj2 = MapP(3:end-1,2:end-1);
    
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/(h/h0)];     % one to the top
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/(h/h0)];     % one to the bottom
    
    
    % Coefficients for x-gradient
    ii  = MapU(2:end-1,2:end-1);
    
    %         left             ||           right
    jj1 = MapP(2:end-1,2:end-2); jj2 = MapP(2:end-1,3:end-1);
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/(h/h0)];     % one to the left
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/(h/h0)];     % one to the right
    
    
    % Assemble coefficient matrix
    GG  = sparse(IIL,JJL,AAL,NW+NU,NP);
end


%% assemble coefficients for divergence operator (DD)

if ~exist('DD','var') || bnchm
    IIL = [];       % equation indeces into A
    JJL = [];       % variable indeces into A
    AAL = [];       % coefficients for A
    
    %internal points
    ii  = MapP(2:end-1,2:end-1);
    
    % coefficients multiplying velocities U, W
    %          left U          ||           right U       ||           top W           ||          bottom W
    jj1 = MapU(2:end-1,1:end-1); jj2 = MapU(2:end-1,2:end); jj3 = MapW(1:end-1,2:end-1); jj4 = MapW(2:end,2:end-1);

    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/(h/h0)];  % U one to the left
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/(h/h0)];  % U one to the right
    IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; aa(:)-1/(h/h0)];  % W one above
    IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; aa(:)+1/(h/h0)];  % W one below

    % Assemble coefficient matrix
    DD  = sparse(IIL,JJL,AAL,NP,NW+NU);
end

%% assemble coefficients for Darcy Flux KF and RF

% qDz  
IIL = [];       % equation indeces into L
JJL = [];       % variable indeces into L
AAL = [];       % coefficients for L
IIR = [];       % equation indeces into R
AAR = [];       % forcing entries for R

% top boundary
ii  = MapW(1,:); jj1 = ii;  
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
aa  = zeros(size(ii));
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii  = MapW(end,:); jj1 = ii; jj2 = MapW(end-1,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
aa  = zeros(size(ii)); 
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% left boundary
ii  = MapW((2:end-1),1); jj1 = ii; jj2 = MapW((2:end-1),2);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapW((2:end-1),end); jj1 = ii; jj2 = MapW((2:end-1),end-1);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];  AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];  AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% internal points for qD_z
ii  = MapW(2:end-1,2:end-1);

% coefficients multiplying Darcy flux
aa  = zeros(size(ii)) + K0./Ksw(2:end-1,:);
IIL = [IIL; ii(:)]; JJL = [JJL; ii(:)];   AAL = [AAL; aa(:)];  % pressure at the centre


% set right hand side in z-direction
rr  = Drhomw(2:end-1,:).*g0./(Drho0*g0);
rr(end,end) = 0;  % avoid conflict with boundary conditions at lower right corner (Div.v = 0)
IIR = [IIR; ii(:)];
AAR = [AAR; rr(:)];

% qDx

% top boundary
ii  = MapU(1,:); jj1 = ii; jj2 = MapU(2,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii  = MapU(end,:); jj1 = ii; jj2 = MapU(end-1,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% left boundary
ii  = MapU((2:end-1),1); jj = ii;
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa  = zeros(size(ii));
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapU((2:end-1),end); jj1 = ii; jj2 = MapU((2:end-1),end-1);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
aa  = zeros(size(ii));
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% internal points
ii  = MapU(2:end-1,2:end-1);

% coefficients multiplying Darcy flux
aa  = zeros(size(ii)) + K0./Ksu(:,2:end-1);
IIL = [IIL; ii(:)]; JJL = [JJL; ii(:)];   AAL = [AAL; aa(:)];  % pressure at the centre


% set right hand side in x-direction
rr  = zeros(size(ii));
IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];


% assemble block matrix and right hand side vector
KF  = sparse(IIL,JJL,AAL,NW+NU,NW+NU);
RF  = sparse(IIR,ones(size(IIR)),AAR);


%% assemble coefficients for Fluid Pressure and right-hand side (KP and RP)

IIL = [];       % equation indeces into L
JJL = [];       % variable indeces into L
AAL = [];       % coefficients for L
IIR = [];       % equation indeces into R
AAR = [];       % forcing entries for R

% Bounday points

% top boundary
ii  = MapP(1,:);  jj1 = ii; jj2 = MapP(2,:);
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii  = MapP(end,:);  jj1 = ii; jj2 = MapP(end-1,:); 
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)]; 

% left boundary  
ii  = MapP((2:end-1),1);  jj1 = ii;  jj2 = MapP((2:end-1),2); 
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapP((2:end-1),end); jj1 = ii;  jj2 = MapP((2:end-1),end-1); 
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% Internal Points
ii  = MapP(2:end-1,2:end-1);
aa  = zeros(size(ii));

% Coefficients multiplying fluid pressure Pf
IIL = [IIL; ii(:)]; JJL = [JJL; ii(:)];   AAL = [AAL; aa(:)];  % pressure at the centre

rr  = VolSrc/(Drho0*u0/h0);
IIR = [IIR; ii(:)];
AAR = [AAR; rr(:)];

KP  = sparse(IIL,JJL,AAL,NP,NP);
RP  = sparse(IIR,ones(size(IIR)),AAR,NP,1);



%% assemble coefficients for compressibility diagonal and right-hand side (KC and RC)
IIL = [];       % equation indeces into L
JJL = [];       % variable indeces into L
AAL = [];       % coefficients for L
IIR = [];       % equation indeces into R
AAR = [];       % forcing entries for R

% Bounday points

% top boundary
ii  = MapP(1,:);  jj1 = ii; jj2 = MapP(2,:);
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii  = MapP(end,:);  jj1 = ii;  jj2 = MapP(end-1,:);
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% left boundary  
ii  = MapP((2:end-1),1);  jj1 = ii;  jj2 = MapP((2:end-1),2);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapP((2:end-1),end); jj1 = ii;  jj2 = MapP((2:end-1),end-1);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];


% Internal Points
ii  = MapP(2:end-1,2:end-1);
aa  = zeros(size(ii)) + e0./zeta;

% Coefficients multiplying compaction pressure p
IIL = [IIL; ii(:)]; JJL = [JJL; ii(:)];   AAL = [AAL; aa(:)];  % pressure at the centre

% RHS
rr = zeros(size(ii));
IIR = [IIR; ii(:)];
AAR = [AAR; rr(:)];

KC = sparse(IIL,JJL,AAL,NP,NP);
RC = sparse(IIR,ones(size(IIR)),AAR,NP,1);

%% assemble coefficients for divergence of matrix mass flux (DD)

IIL = [];       % equation indeces into A
JJL = [];       % variable indeces into A
AAL = [];       % coefficients for A

%internal points
ii  = MapP(2:end-1,2:end-1);

% coefficients multiplying velocities U, W
%          left U          ||           right U       ||           top W           ||          bottom W
jj1 = MapU(2:end-1,1:end-1); jj2 = MapU(2:end-1,2:end); jj3 = MapW(1:end-1,2:end-1); jj4 = MapW(2:end,2:end-1);
rho1 = rhou(:,1:end-1)/Drho0;
rho2 = rhou(:,2:end  )/Drho0;
rho3 = rhow(1:end-1,:)/Drho0;
rho4 = rhow(2:end  ,:)/Drho0;

aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; -rho1(:)/(h/h0)];  % U one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; +rho2(:)/(h/h0)];  % U one to the right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; -rho3(:)/(h/h0)];  % W one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; +rho4(:)/(h/h0)];  % W one below

% Assemble coefficient matrix
DDs  = sparse(IIL,JJL,AAL,NP,NW+NU);

%% assemble coefficients for divergence of relative melt mass flux (DD)

IIL = [];       % equation indeces into A
JJL = [];       % variable indeces into A
AAL = [];       % coefficients for A

%internal points
ii  = MapP(2:end-1,2:end-1);

% coefficients multiplying velocities U, W
%          left U          ||           right U       ||           top W           ||          bottom W
jj1 = MapU(2:end-1,1:end-1); jj2 = MapU(2:end-1,2:end); jj3 = MapW(1:end-1,2:end-1); jj4 = MapW(2:end,2:end-1);
M1  = Mu(:,1:end-1)/Drho0;
M2  = Mu(:,2:end  )/Drho0;
M3  = Mw(1:end-1,:)/Drho0;
M4  = Mw(2:end  ,:)/Drho0;

aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; -M1(:)/(h/h0)];  % U one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; +M2(:)/(h/h0)];  % U one to the right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; -M3(:)/(h/h0)];  % W one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; +M4(:)/(h/h0)];  % W one below

% Assemble coefficient matrix
DDm  = sparse(IIL,JJL,AAL,NP,NW+NU);


%% set pressure fix line
ipx = 2:Nx+1;
ipz = Nz+1;
ip0 = MapP(ipz,ipx);
KP(ip0,:)   = 0;
KP(ip0,ip0) = speye(length(ip0));

%% assemble and scale global coefficient matrix and right-hand side vector

% % Sizes of blocks
% [n1, m1] = size(KV);
% [n2, m2] = size(KF);
% [n3, m3] = size(KC);
% 
% % Total size
% Ntot = n1 + n2 + n3;
% 
% % Preallocate LL as sparse
% if ~exist('total_nnz','var'); total_nnz = nnz(KV) + 2*nnz(GG) + nnz(KF) + nnz(KC) + 2*nnz(DD);  end
% LL = spalloc(Ntot, Ntot, total_nnz);
% 
% % Assign blocks
% LL(1:n1,       1:m1      ) = KV;
% LL(1:n1,    m1+1:m1+m2   ) = GG;
% LL(1:n1, m1+m2+1:m1+m2+m3) = GG;
% 
% LL(n1+1:n1+n2,    1:m1    ) = DD;
% LL(n1+1:n1+n2, m1+1:m1+m2 ) = KF;
% 
% LL(n1+n2+1:n1+n2+n3,       1:m1      ) = DD;
% LL(n1+n2+1:n1+n2+n3, m1+m2+1:m1+m2+m3) = KC;

OV = sparse(NW+NU,NW+NU);
OP = sparse(NW+NU,NP);
OC = sparse(NP,NP);

% Assign blocks
LL = [KV   OV   GG   GG; ...
      OV.' KF   GG   OP; ...
      DDs  DDm  KP   OC; ...
      DD   OP.' OC.' KC];

RR  = [RV; RF; RP; RC];



%% Setting Pc to zero where there is no melt 

bc_ind = [];
bc_ind = [bc_ind;find(twophsw(:)<=0.0) +   NW+NU     ];
bc_ind = [bc_ind;find(twophsu(:)<=0.0) + 2*NW+NU     ];
bc_ind = [bc_ind;find(twophs (:)<=0.0) + 2*NW+2*NU+NP];

bc_val = zeros(size(bc_ind));

% assemble and sort all boundary indices and values
[BC.ind,ind]    =  sort(bc_ind);
 BC.val         =  bc_val(ind);

% extract boundary conditions to reduce problem size
BC.free         =  1:length(LL(1,:));
BC.free(BC.ind) =  [];
TMP             =  LL(:,BC.ind);
RR              =  RR - TMP*BC.val;
LL              =  LL(BC.free,BC.free);
RR              =  RR(BC.free);


%% Scaling coefficient matrix

etagh = ones(Nz+2,Nx+2);  etagh(2:end-1,2:end-1) = eta/e0;
scl = ([zeros(NU+NW,1); zeros(NU+NW,1); 1./etagh(:); zeros(NP,1)]);
scl(BC.ind) = [];
SCL = (abs(diag(LL)) + scl).^0.5;  
SCL = spdiags(1./SCL, 0, length(SCL), length(SCL));

LL = SCL * LL * SCL;
RR = SCL * RR;


%% Solve linear system of equations for W, U, Pf, Pc

if iter==2; pcol = colamd(LL); end               % update column permutation for sparsity pattern once per time step
dLL          = decomposition(LL(:,pcol), 'lu');  % get LU-decomposition for consistent performance of LL \ RR
SS(pcol,1)   = dLL \ RR;                         % solve permuted decomposed system

SOL(BC.free) = SCL * SS;  % update solution
SOL(BC.ind ) = BC.val;    % fill in boundary conditions  
clear SS;                 % clear temporary solution vector

% map solution vector to 2D arrays
W    = full(reshape(SOL(MapW(:))               ,Nz+1,Nx+2));  % matrix z-velocity
U    = full(reshape(SOL(MapU(:))               ,Nz+2,Nx+1));  % matrix x-velocity
wm   = full(reshape(SOL(MapW(:)+(  NW+  NU   )),Nz+1,Nx+2));  % segregation z-velocity
um   = full(reshape(SOL(MapU(:)+(  NW+  NU   )),Nz+2,Nx+1));  % segregation x-velocity
Pf   = full(reshape(SOL(MapP(:)+(2*NW+2*NU   )),Nz+2,Nx+2));  % matrix dynamic pressure
Pc   = full(reshape(SOL(MapP(:)+(2*NW+2*NU+NP)),Nz+2,Nx+2));  % matrix compaction pressure

% redimensionalise solution and parameters
W      = W   *u0;
U      = U   *u0;
wm     = wm *u0;
um     = um *u0;
Pf     = Pf  *p0;
Pc     = Pc  *p0;

if ~bnchm

     % update phase velocities
    Wx  = W + 0.;  % xtl z-velocity
    Ux  = U + 0.;  % xtl x-velocity
    Wm  = W + wm;  % mlt z-velocity
    Um  = U + um;  % mlt x-velocity


end

FMtime = FMtime + toc;
