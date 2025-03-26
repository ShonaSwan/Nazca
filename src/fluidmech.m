tic;

if ~bnchm && step>0 && ~restart

% %***  update mixture mass density (The net compressibility)
% drhodt  = advn_rho;% + (RHO-rho)/dt;
% 
% % residual of mixture mass evolution
% res_rho = (a1*rho-a2*rhoo-a3*rhooo)/dt - (b1*drhodt + b2*drhodto + b3*drhodtoo);
% 
% % volume source and background velocity passed to fluid-mechanics solver
% upd_rho = - res_rho./b1./rho; % + beta*upd_rho;
% VolSrc  = VolSrc + upd_rho;  % correct volume source term by scaled residual

UBG     = - 0*mean(VolSrc,'all')./2 .* (L/2-XXu);
WBG     = - 2*mean(VolSrc,'all')./2 .* (0  -ZZw);

dPchmbdt  = mod_wall*mean(VolSrc,'all') - mod_wall/eta_wall*Pchmb;
res_Pchmb = (a1*Pchmb-a2*Pchmbo-a3*Pchmboo)/dt - (b1*dPchmbdt + b2*dPchmbdto + b3*dPchmbdtoo);

upd_Pchmb = - alpha*res_Pchmb*dt/a1/3 + beta*upd_Pchmb;
Pchmb     = Pchmb + upd_Pchmb;

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
aa  = zeros(size(ii));% + WBG(1,:);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii  = MapW(end,:); jj1 = ii; jj2 = MapW(end-1,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
aa  = zeros(size(ii)) + WBG(end,:);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% left boundary
ii  = MapW(2:end-1,1); jj1 = ii; jj2 = MapW(2:end-1,2);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapW(2:end-1,end); jj1 = ii; jj2 = MapW(2:end-1,end-1);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];  AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];  AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];


% internal points
ii    = MapW(2:end-1,2:end-1);
EtaC1 =  etaco(2:end-1,1:end-1);   EtaC2 =  etaco(2:end-1,2:end);
EtaP1 =  eta  (1:end-1,:      );   EtaP2 =  eta  (2:end,:      );

% coefficients multiplying z-velocities W
%             top          ||         bottom          ||           left            ||          right
jj1 = MapW(1:end-2,2:end-1); jj2 = MapW(3:end,2:end-1); jj3 = MapW(2:end-1,1:end-2); jj4 = MapW(2:end-1,3:end);

aa  = 2/3*(EtaP1+EtaP2)/h^2 + 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % W on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-2/3*EtaP1(:)/h^2];      % W one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-2/3*EtaP2(:)/h^2];      % W one below
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-1/2*EtaC1(:)/h^2];      % W one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-1/2*EtaC2(:)/h^2];      % W one to the right


% what shall we do with the drunken sailor...
 if ~bnchm
     aa  = ddz(rho,h).*g0.*dt;
    IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)];
 end

% coefficients multiplying x-velocities U
%         top left         ||        bottom left          ||       top right       ||       bottom right
jj1 = MapU(2:end-2,1:end-1); jj2 = MapU(3:end-1,1:end-1); jj3 = MapU(2:end-2,2:end); jj4 = MapU(3:end-1,2:end);

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % U one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;+(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % U one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and right


% z-RHS vector
rr  = + (rhofz(2:end-1,:) - mean(rhofz(2:end-1,:),2)) .* g0;
if bnchm; rr = rr + src_W_mms(2:end-1,:); end

IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];


%  assemble coefficients of x-stress divergence

% top boundary
ii  = MapU(1,:); jj1 = ii; jj2 = MapU(2,:);
aa  = zeros(size(ii)) + bnd_spr * 2;
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
ii  = MapU(2:end-1,1); jj = ii;
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa  = zeros(size(ii));
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapU(2:end-1,end); jj1 = ii; jj2 = MapU(2:end-1,end-1);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
aa  = zeros(size(ii));
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% internal points
ii    = MapU(2:end-1,2:end-1);
EtaC1 = etaco(1:end-1,2:end-1);  EtaC2 = etaco(2:end,2:end-1);
EtaP1 = eta  (:      ,1:end-1);  EtaP2 = eta  (:      ,2:end);

% coefficients multiplying x-velocities U
%            left          ||          right          ||           top             ||          bottom
 jj1 = MapU(2:end-1,1:end-2); jj2 = MapU(2:end-1,3:end); jj3 = MapU(1:end-2,2:end-1); jj4 = MapU(3:end,2:end-1);

aa  = 2/3*(EtaP1+EtaP2)/h^2 + 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % U on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-2/3*EtaP1(:)/h^2];      % U one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-2/3*EtaP2(:)/h^2];      % U one to the right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-1/2*EtaC1(:)/h^2];      % U one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-1/2*EtaC2(:)/h^2];      % U one below

% what shall we do with the drunken sailor...
if ~bnchm
    aa  = ddx(rho,h).*g0.*dt;
    IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)];
end

% coefficients multiplying z-velocities W
%         top left         ||        top right          ||       bottom left       ||       bottom right
 jj1 = MapW(1:end-1,2:end-2); jj2 = MapW(1:end-1,3:end-1); jj3 = MapW(2:end,2:end-2); jj4 = MapW(2:end,3:end-1);

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % W one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;+(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % W one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % W one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and right


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
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % one to the top
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % one to the bottom
    
    
    % Coefficients for x-gradient
    ii  = MapU(2:end-1,2:end-1);
    
    %         left             ||           right
    jj1 = MapP(2:end-1,2:end-2); jj2 = MapP(2:end-1,3:end-1);
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % one to the left
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % one to the right
    
    
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
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];  % U one to the left
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];  % U one to the right
    IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; aa(:)-1/h];  % W one above
    IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; aa(:)+1/h];  % W one below

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
ii  = MapP(1,:).';  jj1 = ii; %jj2 = MapP(2,:).';
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
%IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
%topdPdz = -top*(bndmode~=5)*((rhom(1  ,icx)+rhom(2    ,icx))/2-mean(rhofz(1,:)))*g0;
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];% + topdPdz(:)];

% bottom boundary
ii  = MapP(end ,:).';  jj1 = ii;  %jj2 = MapP(end-1,:).'; 
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
%IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
%botdPdz =  bot*(bndmode~=5)*((rhom(end,icx)+rhom(end-1,icx))/2-mean(rhofz(end,:)))*g0;
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)]; % + botdPdz(:)];

% left boundary  
ii  = MapP(2:end-1,1);  jj1 = ii;  jj2 = MapP(2:end-1,2); 
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapP(2:end-1,end  ); jj1 = ii; %jj2 = MapP(2:end-1,end-1);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
%IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% Internal Points
ii  = MapP(2:end-1,2:end-1);
jj1 = MapP(1:end-2,2:end-1); % Above
jj2 = MapP(3:end-0,2:end-1); % Below
jj3 = MapP(2:end-1,1:end-2); % Left
jj4 = MapP(2:end-1,3:end-0); % Right

% Coefficients multiplying darcy flow

KD1 = (KD(icz(1:end-2), :) .* KD(icz(2:end-1), :)).^0.5; % above
KD2 = (KD(icz(2:end-1), :) .* KD(icz(3:end-0), :)).^0.5; % below
KD3 = (KD(: ,icx(1:end-2)) .* KD(:, icx(2:end-1))).^0.5; % left
KD4 = (KD(: ,icx(2:end-1)) .* KD(: ,icx(3:end-0))).^0.5; % right
 
aa = (KD1 + KD2 + KD3 + KD4) / h^2;
IIL = [IIL; ii(:)]; JJL = [JJL; ii(:)];    AAL = [AAL; aa(:)];
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; -KD1(:) / h^2];  % Above
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; -KD2(:) / h^2];  % Below
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; -KD3(:) / h^2];  % Left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; -KD4(:) / h^2];  % Right

% RHS

ii = MapP(2:end-1, 2:end-1);

KD_rho_g     = KD .* (rhom - mean(rho,2)) * g0;
dKD_rho_g_dz = ddz((KD_rho_g([1 1:end],:)+KD_rho_g([1:end end],:))/2,h);
rr  = -dKD_rho_g_dz;
IIR = [IIR; ii(:)];
AAR = [AAR; rr(:)];

KF = sparse(IIL,JJL,AAL,NP,NP);
RF = sparse(IIR,ones(size(IIR)),AAR,NP,1);
% 
% % set P = 0 in fixed point
% nzp = round((Nz+2)/2);
% nxp = round((Nx+2)/2);
% np0 = MapP(nzp,nxp);
%  %KF(np0,np0) = KF(np0,np0)+1e-12;
% KF(np0,:  ) = 0;
% KF(np0,np0) = 1;
% DD(np0,:  ) = 0;
% RF(np0)     = 0;


%% assemble coefficients for compressibility diagonal and right-hand side (KC and RC)
IIL = [];       % equation indeces into L
JJL = [];       % variable indeces into L
AAL = [];       % coefficients for L
IIR = [];       % equation indeces into R
AAR = [];       % forcing entries for R

% Bounday points

% top boundary
ii  = MapP([1,2],:).';  jj1 = ii; %jj2 = MapP(2,:).';
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
%IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii  = MapP(end ,:).';  jj1 = ii;  %jj2 = MapP(end-1,:).'; 
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
%IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% left boundary  
ii  = MapP(2:end-1,1);  jj1 = ii;  jj2 = MapP(2:end-1,2);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapP(2:end-1,end  ); jj1 = ii; %jj2 = MapP(2:end-1,end-1);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
%IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];


% Internal Points
ii  = MapP(2:end-1,2:end-1);
aa  = zeros(size(ii)) + 1./zeta;

% Coefficients multiplying compaction pressure p
IIL = [IIL; ii(:)]; JJL = [JJL; ii(:)];   AAL = [AAL; aa(:)];  % pressure at the centre

% RHS
rr = zeros(size(ii));
IIR = [IIR; ii(:)];
AAR = [AAR; rr(:)];

KC = sparse(IIL,JJL,AAL,NP,NP);
RC = sparse(IIR,ones(size(IIR)),AAR,NP,1);


%% assemble and scale global coefficient matrix and right-hand side vector
OO = zeros(NP, NP);

LL  = [ KV   GG  GG  ; ...
        DD   KF  OO  ; ...
        DD   OO  KC  ];

RR  = [RV; RF; RC];


SCL = (abs(diag(LL))).^0.5;
SCL = diag(sparse( 1./(SCL + sqrt(h^2./geomean(eta(:)))/1000) ));

FF  = LL*[W(:);U(:);Pf(:);Pc(:)] - RR;

LL  = SCL*LL*SCL;
FF  = SCL*FF;
RR  = SCL*RR;


%% Solve linear system of equations for vx, vz, P, Pc

SOL = SCL*(LL\RR);  % update solution

% map solution vector to 2D arrays
W  = full(reshape(SOL(MapW(:))           ,Nz+1,Nx+2));  % matrix z-velocity
U  = full(reshape(SOL(MapU(:))           ,Nz+2,Nx+1));  % matrix x-velocity
Pf = full(reshape(SOL(MapP(:)+(NW+NU   )),Nz+2,Nx+2));  % matrix dynamic pressure
Pc = full(reshape(SOL(MapP(:)+(NW+NU+NP)),Nz+2,Nx+2));  % matrix compaction pressure

% upd_W = - alpha*full(reshape(SOL(MapW(:))        ,Nz+1,Nx+2)) + beta*upd_W;  % matrix z-velocity
% upd_U = - alpha*full(reshape(SOL(MapU(:))        ,Nz+2,Nx+1)) + beta*upd_U;  % matrix x-velocity
% upd_P = - alpha*full(reshape(SOL(MapP(:)+(NW+NU)),Nz+2,Nx+2)) + beta*upd_P;  % matrix dynamic pressure
% 
% W = W + upd_W;
% U = U + upd_U;
% P = P + upd_P;

% end


if ~bnchm

    % z-Darcy flux
    
    qDz(2:end-1,2:end-1) = - (KD(1:end-1,:).*KD(2:end,:)).^0.5.*(ddz(Pf(2:end-1,2:end-1),h)-((rhom(1:end-1,:)+rhom(2:end,:))/2-mean(rhofz(2:end-1,:),2)).*g0); % melt segregation speed
    qDz([1,end],:) = min(1,1-[top;bot]).*qDz([2,end-1],:);
    % qDz([1,end],:) = min(1,1).*qDz([2,end-1],:);
    qDz(:,[1 end]) = qDz(:,[2 end-1]);
 

    qDx(2:end-1,2:end-1) = - (KD(:,1:end-1).*KD(:,2:end)).^0.5 .*(ddx(Pf(2:end-1,2:end-1),h) - ((rhom(:,1:end-1)+rhom(:,2:end))/2 - mean(rhofz(1:end-1,2:end-1), 2)));
    qDx(:,[1,end]) = qDx(:,[2,end-1]); % Simple extrapolation for left/right, adjust if needed
    qDx([1 end],:) = qDx([2 end-1],:); % Top/bottom copied from interior
    

    % % phase segregation speeds
    % wm(2:end-1,2:end-1) = ((rhom(1:end-1,:)+rhom(2:end,:))/2-mean(rhofz(2:end-1,:),2)).*g0.*(Ksgr_m(1:end-1,:).*Ksgr_m(2:end,:)).^0.5; % melt segregation speed
    % wm([1,end],:) = min(1,1-[top;bot]).*wm([2,end-1],:);
    % if periodic
    %     wm(:,[1 end]) = wm(:,[end-1 2]);
    % else
    %     wm(:,[1 end]) = wm(:,[2 end-1]);
    % end
    % 
    % wx(2:end-1,2:end-1) = ((rhox(1:end-1,:)+rhox(2:end,:))/2-mean(rhofz(2:end-1,:),2)).*g0.*(Ksgr_x(1:end-1,:).*Ksgr_x(2:end,:)).^0.5; % solid segregation speed
    % wx([1,end],:) = min(1,1-[top;bot]).*wx([2,end-1],:);
    % if periodic
    %     wx(:,[1 end]) = wx(:,[end-1 2]);
    % else
    %     wx(:,[1 end]) = wx(:,[2 end-1]);
    % end

    chiz = (chi(icz(1:end-1),icx)+chi(icz(2:end),icx))./2;
    chix = (chi(icz,icx(1:end-1))+chi(icz,icx(2:end)))./2;

    muz  = (mu (icz(1:end-1),icx)+mu (icz(2:end),icx))./2;
    mux  = (mu (icz,icx(1:end-1))+mu (icz,icx(2:end)))./2;

    wm = qDz./muz;
    um = qDx./mux;

    % update phase velocities
    % Wf  = W + wf;  % mvp z-velocity
    % Uf  = U + 0.;  % mvp x-velocity
    Wx  = W + 0.;  % xtl z-velocity
    Ux  = U + 0.;  % xtl x-velocity
    Wm  = W + wm;  % mlt z-velocity
    Um  = U + um;  % mlt x-velocity

    
    %% update time step
    dtk = (h/2)^2/max([kc(:);(kT(:)+ks(:).*T(:))./rho(:)./cP(:)]); % diffusive time step size  
    dta =  h/2   /max(abs([Um(:).* mux(:);Wm(:).* muz(:); ...  % advective time step size
                           Ux(:).*chix(:);Wx(:).*chiz(:)])); %+eps
    dtc = maxcmp./max(abs([advn_X(:)./rho(:);advn_M(:)./rho(:)]));
    dt  = min([1.1*dto,min(CFL*[dtk,dta,dtc]),dtmax]);                         % time step size
end

% end

FMtime = FMtime + toc;