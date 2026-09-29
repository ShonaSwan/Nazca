% set directory paths
addpath('./ternplot');

% create output directory
if ~isfolder([outdir,'/',runID])
    mkdir([outdir,'/',runID]);
end

% save input parameters and runtime options (unless restarting)
if restart == 0 && save_op == 1
    parfile = [outdir,'/',runID,'/',runID,'_par'];
    save(parfile);
end

fprintf('\n\n')
fprintf('*************************************************************\n');
fprintf('*****  RUN NAZCA MODEL | %s  **************\n',datetime('now'));
fprintf('*************************************************************\n');
fprintf('\n   run ID: %s \n\n',runID);


%define individualy
                 % load custom colormap
run(['../cal/cal_',calID]);  % load melt model calibration

% load and process custom colormaps
switch colourmap
    case 'ocean'
        load ./colmap/ocean.mat
        colmap = ocean;   
    case 'acton'
        load ./colmap/acton.mat 
        colmap = acton;
    case 'devon'
        load ./colmap/devon.mat
        colmap = devon;
    case 'lajolla'
        load ./colmap/lajolla.mat
        colmap = lajolla;
    case 'lipari'
        load ./colmap/lipari.mat
        colmap = lipari;  
    case 'lapaz'
        load ./colmap/lapaz.mat
        colmap = lapaz;  
    case 'glasgow'
        load ./colmap/glasgow.mat
        colmap = glasgow;  
end
nclmp = length(colmap);
lncls = colmap([5,135],:);


    BCA     =  {'',''};  % boundary condition on advection (top,bot,left,right)
    BCD     =  {'',''};  % boundary condition on diffusion (top,bot,left,right)

Dsm = cal.Dsx;


% normalise major components to anhydrous unit sum, rescale to hydrous
maj0      = maj0     ./(sum(maj0     )+sum(trc0     *trcscale));
maj1      = maj1     ./(sum(maj1     )+sum(trc1     *trcscale));
maj_plume = maj_plume./(sum(maj_plume)+sum(trc_plume*trcscale));
maj_crust = maj_crust./(sum(maj_crust)+sum(trc_crust*trcscale));
dg_maj    = dg_maj-round(mean([dg_maj,dg_trc*trcscale]),16);
dr_maj    = dr_maj-round(mean([dr_maj,dr_trc*trcscale]),16);

% get coordinate arrays
Xc        = -h/2:h:L+h/2;
Zc        = -h/2:h:D+h/2;
Xf        = (Xc(1:end-1)+Xc(2:end))./2; Xu = Xf;
Zf        = (Zc(1:end-1)+Zc(2:end))./2; Zw = Zf;
[XXu,ZZu] = meshgrid(Xf,Zc);
[XXw,ZZw] = meshgrid(Xc,Zf);
Xc        = Xc(2:end-1);
Zc        = Zc(2:end-1);
[XX,ZZ]   = meshgrid(Xc,Zc);

Nx = length(Xc);
Nz = length(Zc);

% get smoothed initialisation field
rng(seed);
smth = smth*Nx*Nz*1e-4;
rp   = randn(Nz,Nx);
for i = 1:round(smth)
    rp = rp + diffus(rp,1/8*ones(size(rp)),1,[1,2],BCD);
    rp = rp - mean(mean(rp));
end
rp = rp./max(abs(rp(:)));

gp = exp(-(XX-L/2).^2/(max(L,D)/8)^2 - (ZZ-D/2).^2/(max(L,D)/8)^2);

% get mapping arrays
NP = (Nz+2) * (Nx+2);
NW = (Nz+1) * (Nx+2);
NU = (Nz+2) * (Nx+1);
NV = NW+NU;
MapP = reshape(1:NP,Nz+2,Nx+2);
MapW = reshape(1:NW,Nz+1,Nx+2);
MapU = reshape(1:NU,Nz+2,Nx+1) + NW;


% set up shape functions for transient boundary layers
topshape = zeros(size(ZZ));
botshape = zeros(size(ZZ));

if ~any(bnd_h)
    switch bndmode
        case 0  % Mid ocean ridge set up
            topshape = exp( ( -ZZ+h/2)/bnd_w);
            botshape = exp(-(D-ZZ-h/2)/bnd_w);
   
        case 1  % Plume 
           topshape = exp( ( -ZZ+h/2)/bnd_w);
           botshape = exp(-(D-ZZ-h/2)/bnd_w);

        case 2  % Plume-Ridge set up
           topshape = exp( ( -ZZ+h/2)/bnd_w);
           botshape = exp(-(D-ZZ-h/2)/bnd_w);
    end
end 

bnd_S = zeros(Nz,Nx);
bnd_C = zeros(Nz,Nx,cal.ncmp+cal.ntrc);

% Boundary condition options:
% -1 = free slip, +1 = No slip, 0 = No gradient/Symmetrical
if bndmode == 0 % Mid Ocean Ridge 1/2 spreading 
Wtop   =  0;  Wbot  = -1;  Wleft  = -1;  Wright  = -1;
Utop   = +1;  Ubot  = -1;  Uleft  =  0;  Uright  = -1;
qDxtop = -1; 
Pall = -1;  
elseif bndmode == 1 % Mid Ocean Ridge full ridge 
Wtop   =  0;  Wbot  = -1;  Wleft  = -1;  Wright  = -1;
Utop   = -1;  Ubot  = -1;  Uleft  = -1;  Uright  = -1;
qDxtop = -1; 
Pall = -1;  
elseif bndmode == 2 % Plume only (closed sides)
Wtop  =  0;   Wbot  = -1;  Wleft  = -1;  Wright  = -1;
Utop  =  -1;  Ubot  = -1;  Uleft  =  0;  Uright  =  0;
qDxtop = -1;
Pall = -1;
else %Error Message 
 disp('Invalid Bndmode');
end

%Se ghosted index array
    icx = [1,1:Nx,Nx];
    icz = [1,1:Nz,Nz];
    ifx = [2,1:Nx+1,Nx];
    ifz = [2,1:Nz+1,Nz];

% initialise solution fields

% Calculating the crustal thickness 
if crust_sw == 0
   Hc = Hcmin - 1e6;
elseif crust_sw == 1
   Hc = Hcmin - max(0,7e3 - Hcmin) * (1 - exp(-200 * (sprate * yr)));
else
   disp('Invalid crust setting'); 
end 

switch init_mode
    case 'plume'
        pl_sigma = pl_width / (2 * sqrt(2 * log(2)));
        Tp = T0 + (T1 - T0) * erf(ZZ ./ (2 * sqrt(1e-6 * minage)));

        %Defining the Gaussian inflow profile of the plume 
        pl_profile = exp(-((XX - pl_local).^2) / pl_width^2 - ((ZZ - D).^2) / (2*pl_width)^2); % Peaks at z=D, x=pl_local

        Tp = Tp + dT_plume .* pl_profile + dTr.*rp + dTg.*gp;
        
        maj = zeros(Nz,Nx,cal.ncmp);
        for i = 1:cal.ncmp
          maj(:,:,i) = maj0(i) + (maj_plume(i) - maj0(i)) .* pl_profile + (maj_crust(i) - maj0(i)) .* (1 - erf((ZZ/D - Hc/D + rp*h*dlay)/wlay_c))/2 + dr_maj(i).*rp + dg_maj(i).*gp;
        end
            
        trc = zeros(Nz,Nx,cal.ntrc);
        for i = 1:cal.ntrc
         trc(:,:,i) = trc0(i) + (trc_plume(i) - trc0(i)) .* pl_profile + (trc_crust(i) - trc0(i)) .* (1 - erf((ZZ/D - Hc/D + rp*h*dlay)/wlay_c))/2 + dr_trc(i).*rp + dg_trc(i).*gp;
        end
    case 'MOR'
        if bndmode == 0
            sprtime = XX./sprate + minage;
        elseif bndmode == 1
            sprtime = abs(XX-L/2)./sprate + minage;
        end
        Tp = T0 + (T1 - T0) * erf(ZZ ./ (2 * sqrt(1e-6 * sprtime)));

        maj = zeros(Nz,Nx,cal.ncmp);
        for i = 1:cal.ncmp
            maj(:,:,i)  =  maj0(i) + (maj_crust(i)-maj0(i)) .* (1-erf((ZZ/D-Hc/D+rp*h*dlay)/wlay_c))/2 + dr_maj(i).*rp + dg_maj(i).*gp;  % major elements
        end

        trc = zeros(Nz,Nx,cal.ntrc);
        for i = 1:cal.ntrc
            trc(:,:,i)  =  trc0(i) + (trc_crust(i)-trc0(i)) .* (1-erf((ZZ/D-Hc/D+rp*h*dlay)/wlay_c))/2 + dr_trc(i).*rp + dg_trc(i).*gp;  % trace elements
        end

    case 'PRI'
        if bndmode == 0
            sprtime = XX./sprate + minage;
        else
            sprtime = abs(XX-L/2)./sprate + minage;
        end
        Tp_MOR = T0 + (T1 - T0) * erf(ZZ ./ (2 * sqrt(1e-6 * sprtime)));

        pl_sigma = pl_width / (2 * sqrt(2 * log(2)));
        pl_profile = exp(-((XX - pl_local).^2) / pl_width^2 - ((ZZ - D).^2) / (2*pl_width)^2); % Peaks at z=D, x=pl_local
        Tp_plume = dT_plume .* pl_profile;
        
        Tp = Tp_MOR + Tp_plume + dTr.*rp + dTg.*gp;
              
        maj = zeros(Nz,Nx,cal.ncmp);
        for i = 1:cal.ncmp
            maj(:,:,i)  =  maj0(i) + ...
            (maj_plume(i) - maj0(i)) .* pl_profile + ...
            (maj_crust(i) - maj0(i)) .* (1-erf((ZZ/D-Hc/D+rp*h*dlay)/wlay_c))/2 + dr_maj(i).*rp + dg_maj(i).*gp;  % major elements
        end

        trc = zeros(Nz,Nx,cal.ntrc);
        for i = 1:cal.ntrc
            trc(:,:,i) = trc0(i) + ...
            (trc_plume(i) - trc0(i)) .* pl_profile + ...
            (trc_crust(i) - trc0(i)) .* (1-erf((ZZ/D-Hc/D+rp*h*dlay)/wlay_c))/2 + dr_trc(i).*rp + dg_trc(i).*gp;  % trace elements
        end
    
end 

% combine major and trace components
c = cat(3,maj,trc*trcscale); c = c./sum(c,3);

%Defining the top bounday spreading rate 's' shape function
%bnd_spr = (1-exp(-Xu./bnd_sprw)) .* sprate;
% apply initial boundary layers

Tin   = Tp;
majin = maj;
trcin = trc;

U   =  zeros(Nz+2,Nx+1);  UBG = U; Ui = U; upd_U = 0*U;  um = 0.*U; qDx = 0.*U; Umix = 0.*U;
W   =  zeros(Nz+1,Nx+2);  WBG = W; Wi = W; wx = 0.*W; wm = 0.*W; upd_W = 0*W;  qDz = 0.*W; Wmix = 0.*W;
Pf  =  zeros(Nz+2,Nx+2);  Vx = 0.*Tp; Vm = 0.*Tp; upd_Pf= 0*Pf; %Div_rhoV = 0.*P;  DD = sparse(length(P(:)),length([W(:);U(:)]));
Pc  =  zeros(Nz+2,Nx+2);  P = Pf + Pc;
SOL = [W(:);U(:);wm(:);um(:);Pf(:);Pc(:)]; 

% initialise auxiliary fields
Wx  = W;  Ux  = U;
Wm  = W;  Um  = U;

% define dimensional scales
h0    = D/10;
e0    = 1e18;
Drho0 = 500;
p0    = Drho0*g0*h0;
u0    = p0*h0./e0;
K0    = h0^2/e0;
t0    = h0/u0;
MFD0  = Drho0/t0;
CMP0  = 1/t0;

% initialise other parameters
Re     = eps;
Div_V  = 0.*Tp;  advn_rho = 0.*Tp;  advn_X = 0.*Tp; advn_M = 0.*Tp;
drhodt = 0.*Tp;  drhodto = drhodt; 
exx    = 0.*Tp;  ezz = 0.*Tp;  exz = zeros(Nz-1,Nx-1);  eII = 0.*Tp;
txx    = 0.*Tp;  tzz = 0.*Tp;  txz = zeros(Nz-1,Nx-1);  tII = 0.*Tp;
eta    = 1e21.*ones(Nz,Nx);  etai = eta;
zeta   = 100.*eta;  zetai = zeta;
etamax = etacntr.*etamin;
KD     = 1e-24.*ones(Nz,Nx); KDi = KD;
MFDSrc = 0.*Tp;
CMPSrc = 0.*Tp;
MFDCrr = 0.*Tp;
CMPCrr = 0.*Tp;
kW     = 0.*Tp;
Tref   = min(cal.T0) + 273.15;
Pref   = 1e5;
sref   = 0e3; % reference entropy 
rhom0  = mean(cal.rhox0-300).*ones(size(Tp));
rhox0  = mean(cal.rhox0).*ones(size(Tp)); 
Pt     = Ptop + mean(rhox0,'all').*g0.*ZZ;  Pl = Pt;  Pto = Pt; Ptoo = Pt; dPtdt = 0*Pt; dPtdto = dPtdt; dPtdtoo = dPtdto;
rhox   = rhox0.*(1-aTx.*(Tp-Tref)+bPx.*(Pt-Pref));
rhom   = rhom0.*(1-aTx.*(Tp-Tref)+bPm.*(Pt-Pref));
rho    = rhox;
M      = 0*rhox;
X      = rhox;
rhow   = (rho(icz(1:end-1),:)+rho(icz(2:end),:))/2;
rhou   = (rho(:,icx(1:end-1))+rho(:,icx(2:end)))/2;
rhoWo  = rhow.*W(:,2:end-1); rhoWoo = rhoWo; advn_mz = 0.*rhoWo(2:end-1,:);
rhoUo  = rhou.*U(2:end-1,:); rhoUoo = rhoUo; advn_mx = 0.*rhoUo;
deltac = ones(Nz,Nx);
mq     = zeros(size(Tp));  xq = 1-mq;
ups    = zeros(size(Tp));
cal.Ktrc = ones(Nz*Nx,cal.ntrc);

aT   = aTx;
kT   = kTx;
cP   = cPx; RhoCp = rho.*cP;
Adbt = aT./RhoCp;
Tp   = (Tp+273.15); %T = Tp;
T    = Tp.*exp(Adbt.*(Pt-Pref));
diss = 0.*T;
x    = xq;  m = mq; mu = m; chi = x;
dto  = dt;
step = 0;

% get volume fractions and bulk density
%step    = 0;
EQtime  = 0;
FMtime  = 0;
TCtime  = 0;
UDtime  = 0;
a1      = 1; a2 = 1; a3 = 0; b1 = 1; b2 = 0; b3 = 0;

res  = 1;  tol = 1e-7;  it = 1; iter = 1;

while res > tol

    Pti = Pt; Ti = T; xi = xq;          

    % %%%%The melt model %%%%

        eqtime = tic;

        var.c      = reshape(c,Nx*Nz,cal.ncmp+cal.ntrc);  % component fractions [wt]
        var.T      = reshape(T,Nx*Nz,1)-273.15;   % temperature [C]
        var.P      = reshape(Pt,Nx*Nz,1)/1e9;     % pressure [GPa]
        var.m      = reshape(mq,Nx*Nz,1);         % melt fraction [wt](melt model 

        [var,cal] = leappart(var,cal,'E');

        Tsol   = reshape(cal.Tsol,Nz,Nx);
        Tliq   = reshape(cal.Tliq,Nz,Nx);

        mq = reshape(var.m.*(var.m>eps^0.5),Nz,Nx);
        xq = reshape(var.x.*(var.x>eps^0.5),Nz,Nx);  
        mq = mq./(mq+xq);
        xq = xq./(mq+xq);
        x  = xq;  m = mq;

        cxq  = reshape(var.cx,Nz,Nx,cal.ncmp+cal.ntrc); cxq = cxq./sum(cxq,3);
        cmq  = reshape(var.cm,Nz,Nx,cal.ncmp+cal.ntrc); cmq = cmq./sum(cmq,3);
        cm   = cmq;
        cx   = cxq;
        c    = m.*cm + x.*cx;        c    = c   ./sum(c,3);
        majm = cmq(:,:,1:cal.ncmp);  majm = majm./sum(majm,3);
        majx = cxq(:,:,1:cal.ncmp);  majx = majx./sum(majx,3);
        maj  = m.*majm + x.*majx;    maj  = maj ./sum(maj ,3);
        trcm = cmq(:,:,cal.ncmp+1:end)/trcscale; 
        trcx = cxq(:,:,cal.ncmp+1:end)/trcscale; 
        trc  = m.*trcm + x.*trcx;

        update;

        Ktrc  = zeros(Nz,Nx,cal.ntrc);
        for i = 1:cal.ntrc
            % update bulk partitioning coefficients
            for j=1:cal.nmem; Ktrc(:,:,i) = Ktrc(:,:,i) + cal.Ktrc_mem(i,j) .* majx_mem(:,:,j)./100; end
        end
        cal.Ktrc = reshape(Ktrc,Nz*Nx,cal.ntrc);

        sref = 0;
        sm   = cPm.*log(Tp./Tref) + Dsm;  
        sx   = cPx.*log(Tp./Tref);


        % Removing melt to get a suitable initial melt fraction
        if it>10 && any(m(:)>minit)
            mi  = m;
            if bndmode == 0
                rdist = max(1e-6,min(1-1e-6,sqrt((XX./L).^2 + (ZZ./D).^2)));
            else
                rdist = max(1e-6,min(1-1e-6,sqrt((XX./L*2-1).^2 + (ZZ./D).^2)));
                for i=1:5; rdist = rdist + diffus(rdist,1/8*ones(size(rdist)),1,[1,2],BCD); end
            end
            msc = (minit.^-4 + (m+eps).^-4).^(-1/4).* (1 - rdist).^(1/2);
            m   = m .* (msc./(m+eps)).^0.01;
            % m(end-1:end,:) = 0;
            % m   = m + diffus(m,1e-3*ones(size(rp)),1,[1,2],BCD);

            SUM = x+m;
            x = x./SUM;  m = m./SUM;

            c    = x.*cxq + m.*cmq;
            Tp   = Tp  - (mi-m).*T.*Dsm./cP;
            sm   = cPm.*log(Tp./Tref) + Dsm;
            sx   = cPx.*log(Tp./Tref);
            s    = x.*sx + m.*sm;
            rho  = 1./(m./rhom  + x./rhox);
        else
            s = x.*sx + m.*sm;
        end

        X    = rho.*x;
        M    = rho.*m;  RHO = X+M;
        Cm   = M.*cm;
        Cx   = X.*cx;
        C    = Cm + Cx;
        Sm   = M.*sm;
        Sx   = X.*sx;
        S    = Sm + Sx;

        [T ,si] = StoT(T ,S./rho,cat(3,Pt,Pt),cat(3,m,x),[cPm;cPx],[aTm;aTx],[bPm;bPx],cat(3,rhom0,rhox0),[sref+Dsm;sref],Tref,Pref);

        res  = norm(Pt(:)-Pti(:),2)./norm(Pt(:),2) ...
             + norm( T(:)-Ti  (:),2)./norm( T(:),2);

        it = it+1;
end

Pto  = Pt;
To   = T;  
Tpo  = Tp;
So   = S;  
Mo   = M;
Xo   = X;
Cmo  = Cm;
Cxo  = Cx;
rhoo = rho;
rhoxo= rhox;
chio = chi;
muo  = mu;
clear Twall majwall trcwall;
Twall(1,:) = T(1  ,:)-273.15;
Twall(2,:) = T(end,:)-273.15;
majwall(1,:,:) = nan*maj(1,:,:);
majwall(2,:,:) = maj(end,:,:);
trcwall(1,:,:) = nan*trc(1,:,:);
trcwall(2,:,:) = trc(end,:,:)*trcscale;
cwall          = cat(3,majwall,trcwall);

% initialise phase change rates
Gx   = 0.*x; Gm  = 0.*m; 
Gem  = 0.*m; Gex = 0.*x; Gin = 0.*x;
Gemc = 0.*maj; Gexc = 0.*maj;
Gemt = 0.*trc; Gext = 0.*trc;

% initialise auxiliary variables
res_MFD = 0.*rho;  res_CMP = 0.*rho;
dSdt   = 0.*T;  dSdto  = dSdt; diss_h = 0.*T;  res_S = 0.*S;
dCmdt  = 0.*c;  dCmdto  = dCmdt; res_Cm = 0.*cm;
dCxdt  = 0.*c;  dCxdto  = dCxdt; res_Cx = 0.*cx;
qz_advn_Sx = 0.*W; qx_advn_Sx = 0.*U;
qz_advn_Sm = 0.*W; qx_advn_Sm = 0.*U;
qz_diff_S  = 0.*W; qx_diff_S  = 0.*U;
qz_advn_Cx = 0.*W; qx_advn_Cx = 0.*U;
qz_advn_Cm = 0.*W; qx_advn_Cm = 0.*U;
qz_diff_Cx = 0.*W; qx_diff_Cx = 0.*U;
qz_diff_Cm = 0.*W; qx_diff_Cm = 0.*U;

specrad.S.est   = 0.9.*ones(Nz*Nx,1);          specrad.S.mean   = 0.9;
specrad.MFD.est = 0.9.*ones(Nz*Nx,1);          specrad.MFD.mean = 0.9;
specrad.CMP.est = 0.9.*ones(Nz*Nx,1);          specrad.CMP.mean = 0.9;
specrad.Cm.est  = 0.9.*ones(Nz*Nx*(cal.ncmp+cal.ntrc),1); specrad.Cm.mean  = 0.9;
specrad.Cx.est  = 0.9.*ones(Nz*Nx*(cal.ncmp+cal.ntrc),1); specrad.Cx.mean  = 0.9;

GHST.S   = zeros(Nz*Nx, itpar.aa.m+1);
GHST.MFD = zeros(Nz*Nx, itpar.aa.m+1);
GHST.CMP = zeros(Nz*Nx, itpar.aa.m+1);
GHST.Cm  = zeros(Nz*Nx*(cal.ncmp+cal.ntrc), itpar.aa.m+1);
GHST.Cx  = zeros(Nz*Nx*(cal.ncmp+cal.ntrc), itpar.aa.m+1);

FHST.S   = zeros(Nz*Nx, itpar.aa.m+1);
FHST.MFD = zeros(Nz*Nx, itpar.aa.m+1);
FHST.CMP = zeros(Nz*Nx, itpar.aa.m+1);
FHST.Cm  = zeros(Nz*Nx*(cal.ncmp+cal.ntrc), itpar.aa.m+1);
FHST.Cx  = zeros(Nz*Nx*(cal.ncmp+cal.ntrc), itpar.aa.m+1);

% initialise timing and iterative parameters
frst    = 1;
step    = 0;
time    = 0;
hist    = [];
dsumSdto = 0;  dsumSdt = 0;
dsumBdto = 0;  dsumBdt = 0;
dsumMdto = 0;  dsumMdt = 0;
dsumXdto = 0;  dsumXdt = 0;
dsumCmdto = 0;  dsumCmdt = 0;
dsumCxdto = 0;  dsumCxdt = 0;
dsumCdto  = 0;  dsumCdt  = 0;

% tracer switch 
if tracer_sw == 1
num_tracers = 15;
tracer_X = linspace(Xc(1), Xc(end), num_tracers);  % evenly spaced along X
tracer_Z = D * ones(1, num_tracers);               % at the bottom (Z = D)

u_tracer = interp2(XXu, ZZu, U, tracer_X, tracer_Z, 'linear');  
w_tracer = interp2(XXw, ZZw, W, tracer_X, tracer_Z, 'linear');  

tracer_X = tracer_X + u_tracer * dt;
tracer_Z = tracer_Z + w_tracer * dt;
end

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = [outdir,'/',runID,'/',runID,'_cont.mat'];
    elseif restart > 0  % restart from specified continuation frame
        name = [outdir,'/',runID,'/',runID,'_',num2str(restart),'.mat'];
    end
    
    if exist(name,'file')
        fprintf('\n   restart from %s \n\n',name);
        load(name,'U','W','Pf','Pc','Pt','x','m','chi','mu','X','M','S','C','Cm','Cx','T','Tp','c','cm','cx','maj','majm','majx','trc','trcm','trcx','dSdt','dCmdt','dCxdt','drhodt','Gx','Gm','Gem','Gex','Gin','rho','rhox','eta','zeta','Ks','kd','ups','eII','tII','dt','time','step','MFDSrc','MFDCrr','CMPSrc','CMPCrr','wx','wm','cal','specrad');
        name = [outdir,'/',runID,'/',runID,'_hist'];
        load(name,'hist');

        SOL = [W(:);U(:);wm(:);um(:);Pf(:);Pc(:)];
        RHO = X+M; 
     
        store;
        phseql;
        update;
        fluidmech;
        update;
        output;

        time    = time+dt;
        step    = step+1;

    else % continuation file does not exist, start from scratch
        fprintf('\n   !!! restart file does not exist !!! \n   => starting run from scratch %s \n\n',runID);
        store;
        phseql;
        fluidmech;
        update;
        fluidmech;
        update;
        history;
        output;
    end
else
    % complete, plot, and save initial condition
    store;
    phseql;
    fluidmech;
    update;
    fluidmech;
    update;
    history;
    output;
    step = step+1;
end

restart = 0;
