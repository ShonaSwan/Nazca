tic;

% define dimensional scales
h0    = D/10;
e0    = 1e18;
Drho0 = 500;
p0    = Drho0*g0*h0;
u0    = Drho0*g0*h0^2./e0;
KD0   = h0^2/e0;
t0    = h0/u0;

if ~bnchm && step>0 && ~restart

% residual of mixture mass evolution
res_rho = (a1*rho-a2*rhoo-a3*rhooo)/dt - (b1*drhodt + b2*drhodto + b3*drhodtoo);

% volume source and background velocity passed to fluid-mechanics solver
upd_DV  = - alpha*res_rho./b1./rho;
VolSrc  = VolSrc + upd_DV;  % correct volume source term by scaled residual

end


%% assemble coefficients for matrix velocity diagonal and right-hand side (KV and RV)

IIL = [];       % equation indeces into L
JJL = [];       % variable indeces into L
AAL = [];       % coefficients for L
IIR = [];       % equation indeces into R
AAR = [];       % forcing entries for R


% assemble coefficients of z-stress divergence

% top boundary
ii  = MapW(1,:); jj1 = ii; jj2 = MapW(2,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Wtop];
aa  = zeros(size(ii));
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
if bndmode == 1 % Plume set up 
    ii  = MapW(end,:); jj1 = ii; jj2 = MapW(end-1,:);
    aa  = zeros(size(ii)) + bnd_pl(icx) / pl_rate * 2;
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Wbot];
    %aa  = zeros(size(ii)); 
    IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];
else 
    ii  = MapW(end,:); jj1 = ii; jj2 = MapW(end-1,:);
    aa = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Wbot];
    aa  = zeros(size(ii)); 
    IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];
end 

% left boundary
ii  = MapW((2:end-1),1); jj1 = ii; jj2 = MapW((2:end-1),2);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Wleft];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapW((2:end-1),end); jj1 = ii; jj2 = MapW((2:end-1),end-1);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];  AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];  AAL = [AAL; aa(:)+Wright];
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
rr  = + ((rhow(2:end-1,:)-mean(rhow(2:end-1,:),2)).*g0)/(Drho0*g0);
if bnchm; rr = rr + src_W_mms(2:end-1,:); end

IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];


%  assemble coefficients of x-stress divergence

% top boundary
if bndmode == 0 % Mid ocean Ridge set up 
    ii  = MapU(1,:); jj1 = ii; jj2 = MapU(2,:);
    aa  = zeros(size(ii)) + bnd_spr/u0 * 2;
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Utop];
    IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];
else
    ii  = MapU(1,:); jj1 = ii; jj2 = MapU(2,:);
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Utop];
    IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];    
end

% bottom boundary
ii  = MapU(end,:); jj1 = ii; jj2 = MapU(end-1,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Ubot];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% left boundary
ii  = MapU((2:end-1),1); jj = ii; jj2 = MapU((2:end-1),2);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Uleft];
aa  = zeros(size(ii));% + UBG((2:end-1),1);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapU((2:end-1),end); jj1 = ii; jj2 = MapU((2:end-1),end-1);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Uright];
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

%% assemble coefficients for matrix Darcy flow diagonal and right-hand side (KF and RF)

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
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Pftop];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii  = MapP(end,:);  jj1 = ii; jj2 = MapP(end-1,:); 
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Pfbot];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)]; 

% left boundary  
ii  = MapP((2:end-1),1);  jj1 = ii;  jj2 = MapP((2:end-1),2); 
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Pfleft];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapP((2:end-1),end); jj1 = ii;  jj2 = MapP((2:end-1),end-1); 
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Pfright];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% Internal Points
ii  = MapP(2:end-1,2:end-1);
jj1 = MapP(1:end-2,2:end-1); % Above
jj2 = MapP(3:end-0,2:end-1); % Below
jj3 = MapP(2:end-1,1:end-2); % Left
jj4 = MapP(2:end-1,3:end-0); % Right

% Coefficients multiplying darcy flow
KD1 = (KD(icz(1:end-2), :) .* KD(icz(2:end-1), :)).^0.5/KD0; % above
KD2 = (KD(icz(2:end-1), :) .* KD(icz(3:end-0), :)).^0.5/KD0; % below
KD3 = (KD(: ,icx(1:end-2)) .* KD(:, icx(2:end-1))).^0.5/KD0; % left
KD4 = (KD(: ,icx(2:end-1)) .* KD(: ,icx(3:end-0))).^0.5/KD0; % right
 
aa = (KD1 + KD2 + KD3 + KD4)/(h/h0)^2;
IIL = [IIL; ii(:)]; JJL = [JJL; ii(:)];    AAL = [AAL; aa(:)];
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; -KD1(:)/(h/h0)^2];  % Above
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; -KD2(:)/(h/h0)^2];  % Below
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; -KD3(:)/(h/h0)^2];  % Left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; -KD4(:)/(h/h0)^2];  % Right

% RHS

ii = MapP(2:end-1, 2:end-1);

KD_rho_g     = KD .* (rhox - rhom) * g0/(KD0*Drho0*g0);
dKD_rho_g_dz = ddz((KD_rho_g([1 1:end],:)+KD_rho_g([1:end end],:))/2,h/h0);
rr  = dKD_rho_g_dz + VolSrc*t0;
if bndmode == 0
rr(end,end) = 0;  % avoid conflict with boundary conditions at lower right corner (Div.v = 0)
end
IIR = [IIR; ii(:)];
AAR = [AAR; rr(:)];



KF = sparse(IIL,JJL,AAL,NP,NP);
RF = sparse(IIR,ones(size(IIR)),AAR,NP,1);


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
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Pcall];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii  = MapP(end,:);  jj1 = ii;  jj2 = MapP(end-1,:);
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Pcall];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% left boundary  
ii  = MapP((2:end-1),1);  jj1 = ii;  jj2 = MapP((2:end-1),2);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Pcall];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapP((2:end-1),end); jj1 = ii;  jj2 = MapP((2:end-1),end-1);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+Pcall];
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


%% assemble and scale global coefficient matrix and right-hand side vector

% Sizes of blocks
[n1, m1] = size(KV);
[n2, m2] = size(KF);
[n3, m3] = size(KC);

% Total size
Ntot = n1 + n2 + n3;

% Preallocate LL as sparse
if ~exist('total_nnz','var'); total_nnz = nnz(KV) + 2*nnz(GG) + nnz(KF) + nnz(KC) + 2*nnz(DD);  end
LL = spalloc(Ntot, Ntot, total_nnz);

% Assign blocks
LL(1:n1,       1:m1      ) = KV;
LL(1:n1,    m1+1:m1+m2   ) = GG;
LL(1:n1, m1+m2+1:m1+m2+m3) = GG;

LL(n1+1:n1+n2,    1:m1    ) = DD;
LL(n1+1:n1+n2, m1+1:m1+m2 ) = KF;

LL(n1+n2+1:n1+n2+n3,       1:m1      ) = DD;
LL(n1+n2+1:n1+n2+n3, m1+m2+1:m1+m2+m3) = KC;

RR  = [RV; RF; RC];



%% Setting Pc to zero where there is no melt 

bc_ind = find(twophs(:)<=0.0) + NW+NU+NP;
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
scl = ([zeros(NU+NW,1); 1./etagh(:); zeros(NP,1)]);
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
W  = full(reshape(SOL(MapW(:))           ,Nz+1,Nx+2));  % matrix z-velocity
U  = full(reshape(SOL(MapU(:))           ,Nz+2,Nx+1));  % matrix x-velocity
Pf = full(reshape(SOL(MapP(:)+(NW+NU   )),Nz+2,Nx+2));  % matrix dynamic pressure
Pc = full(reshape(SOL(MapP(:)+(NW+NU+NP)),Nz+2,Nx+2));  % matrix compaction pressure

% redimensionalise solution and parameters
W      = W *u0;
U      = U *u0;
Pf     = Pf*p0;
Pc     = Pc*p0;


if ~bnchm

    % z-Darcy flux
    qDz(2:end-1,2:end-1) = - (KD(1:end-1,:).*KD(2:end,:)).^0.5.*(ddz(Pf(2:end-1,2:end-1),h)-((rhom(1:end-1,:)+rhom(2:end,:))/2-mean(rhow(2:end-1,:),2)).*g0); % melt segregation speed
    qDz([1,end],:) = min(1,1-[1;-1]).*qDz([2,end-1],:);
    qDz(:,[1 end]) = qDz(:,[2 end-1]);

    % x-Darcy flux
    qDx(2:end-1,2:end-1) = - (KD(:,1:end-1).*KD(:,2:end)).^0.5 .*(ddx(Pf(2:end-1,2:end-1),h));
    qDx(:,end) = qDx(:,end-1); 
    qDx(:,1) = 0;
    qDx([1 end],:) = qDx([2 end-1],:); % Top/bottom copied from interior

    muz  = (mu (icz(1:end-1),icx)+mu (icz(2:end),icx))./2;
    mux  = (mu (icz,icx(1:end-1))+mu (icz,icx(2:end)))./2;

    wm = qDz./max(mulim,muz);
    um = qDx./max(mulim,mux);

    % update phase velocities
    Wx  = W + 0.;  % xtl z-velocity
    Ux  = U + 0.;  % xtl x-velocity
    Wm  = W + wm;  % mlt z-velocity
    Um  = U + um;  % mlt x-velocity

end

FMtime = FMtime + toc;
