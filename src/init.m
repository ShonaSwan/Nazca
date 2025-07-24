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
load ocean;                  % load custom colormap
run(['../cal/cal_',calID]);  % load melt model calibration


    BCA     =  {'',''};  % boundary condition on advection (top,bot,left,right)
    BCD     =  {'',''};  % boundary condition on diffusion (top,bot,left,right)

Dsm = cal.Dsx;

% normalise major components to anhydrous unit sum, rescale to hydrous
c0(1:end-1) = c0(1:end-1)./sum(c0(1:end-1)).*(1-c0(end));
c1(1:end-1) = c1(1:end-1)./sum(c1(1:end-1)).*(1-c1(end));
cwall(:,1:end-1) = cwall(:,1:end-1)./sum(cwall(:,1:end-1),2).*(1-cwall(:,end));
dcg   = dcg-round(mean(dcg),16);
dcr   = dcr-round(mean(dcr),16);

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
            botshape = exp(-(D-ZZ+h/2)/bnd_w);
   
        case 1  % top only
           topshape = exp( ( -ZZ)/bnd_w);
           botshape = exp(-(D-ZZ)/bnd_w);
    end
end 

bnd_S = zeros(Nz,Nx);
bnd_C = zeros(Nz,Nx,cal.ncmp);
bnd_V = zeros(Nz,Nx);

% Boundary condition options:
% -1 = free slip, +1 = No slip, 0 = No gradient/Symetrical
if bndmode == 0; %Mid Ocean Ridge
Wtop   =  0;  Wbot  = -1;  Wleft  = -1;  Wright  = -1;
Utop   = +1;  Ubot  = -1;  Uleft  =  0;  Uright  = -1;
qDxtop = -1; 
Pall = -1;  

elseif bndmode == 1; % Mantle Plume
Wtop  =  0;   Wbot  = -1;  Wleft  = -1;  Wright  = -1;
Utop  =  -1;  Ubot  = -1;  Uleft  =  0;  Uright  =  0;
qDxtop = -1;
Pall = -1;
% Pftop = -1;  Pfbot = +1;  Pfleft = -1;  Pfright = -1;
% Pcall = -1;

else; %Error Message 
 disp(['Invalid Bndmode']);
end

%Se ghosted index array
    icx = [1,1:Nx,Nx];
    icz = [1,1:Nz,Nz];
    ifx = [2,1:Nx+1,Nx];
    ifz = [2,1:Nz+1,Nz];

% initialise solution fields

% Calculating the crustal thickness 
%Hc = Hcmin - max(0,7e3 - Hcmin) * (1 - exp(-200 * (sprate * yr)));
Hc = Hcmin - 1e6;%max(0,7e3 - Hcmin) * (1 - exp(-200 * (sprate * yr)));   

switch init_mode
    case 'plume'
        pl_sigma = pl_width / (2 * sqrt(2 * log(2)));
        Tp = T0 + (T1 - T0) * erf(ZZ ./ (2 * sqrt(1e-6 * minage)));

        %Defining the Gaussian inflow profile of the plume 
        pl_profile = exp(-((XX - pl_local).^2) / pl_width^2 - ((ZZ - D).^2) / pl_width^2); % Peaks at z=D, x=pl_local

        Tp = Tp + dT_plume .* pl_profile + dTr.*rp + dTg.*gp;
        
        for i = 1:cal.ncmp
            c(:,:,i) = c0(i) + (c1(i) - c0(i)) .* (ZZ/D) .* (1 - pl_profile) + (c_plume(i) - c0(i)) .* pl_profile + dcr(i).*rp + dcg(i).*gp;
        end
        for i = 1:cal.ntrc
            trc(:,:,i) = trc0(i) + (trc1(i) - trc0(i)) .* (ZZ/D) .* (1 - pl_profile) + (trc_plume(i) - trc0(i)) .* pl_profile + dr_trc(i).*rp + dg_trc(i).*gp;
        end

    case 'MOR'
        sprtime = XX./sprate + minage;
        
        Tp = T0 + (T1 - T0) * erf(ZZ ./ (2 * sqrt(1e-6 * sprtime)));

        for i = 1:cal.ncmp
            c(:,:,i)  =  c0(i) + (c_crust(i)-c0(i)) .* (1-erf((ZZ/D-Hc/D+rp*h*dlay)/wlay_c))/2 + dcr(i).*rp + dcg(i).*gp;  % major elements
        end
        trc = zeros(Nz,Nx,cal.ntrc);
        for i = 1:cal.ntrc
            trc(:,:,i)  =  trc0(i) + (trc_crust(i)-trc0(i)) .* (1-erf((ZZ/D-Hc/D+rp*h*dlay)/wlay_c))/2 + dr_trc(i).*rp + dg_trc(i).*gp;  % trace elements
        end

end 


%Defining the top bounday spreading rate 's' shape function
bnd_spr = (1-exp(-Xu./bnd_sprw)) .* sprate;

% apply initial boundary layers

Tin = Tp;
cin = c;
tein = trc;

U   =  zeros(Nz+2,Nx+1);  UBG = U; Ui = U; upd_U = 0*U;  um = 0.*U; qDx = 0.*U;
W   =  zeros(Nz+1,Nx+2);  WBG = W; Wi = W; wx = 0.*W; wm = 0.*W; upd_W = 0*W;  qDz = 0.*W; 
Pf  =  zeros(Nz+2,Nx+2);  Vel = 0.*Tp; upd_Pf= 0*Pf; %Div_rhoV = 0.*P;  DD = sparse(length(P(:)),length([W(:);U(:)]));
Pc   =  zeros(Nz+2,Nx+2);
SOL = [W(:);U(:);Pf(:);Pc(:)];

% initialise auxiliary fields
Wx  = W;  Ux  = U;
Wm  = W;  Um  = U;

Re     = eps;
Div_V  = 0.*Tp;  advn_rho = 0.*Tp;  advn_X = 0.*Tp; advn_M = 0.*Tp;
drhodt = 0.*Tp;  drhodto = drhodt; 
exx    = 0.*Tp;  ezz = 0.*Tp;  exz = zeros(Nz-1,Nx-1);  eII = 0.*Tp;
txx    = 0.*Tp;  tzz = 0.*Tp;  txz = zeros(Nz-1,Nx-1);  tII = 0.*Tp;
eta    = 1e21.*ones(Nz,Nx);
zeta   = 100.*eta;
VolSrc = 0.*Tp;
kW     = 0.*Tp;
Tref   = min(cal.T0) + 273.15;
Pref   = 1e5;
sref   = 0e3; % reference entropy 
c0_oxd = c0*cal.cmp_oxd;
c0_oxd_all = zeros(size(c0,1),9);
c0_oxd_all(:,cal.ioxd) = c0_oxd;
rhom0  = mean(cal.rhox0-500).*ones(size(Tp));
rhox0  = mean(cal.rhox0).*ones(size(Tp)); 
Pchmb  = Pchmb0;  Pchmbo = Pchmb;  Pchmboo = Pchmbo;  dPchmbdt = Pchmb;  dPchmbdto = dPchmbdt; dPchmbdtoo = dPchmbdto;  upd_Pchmb = dPchmbdt;
Pt     = Ptop + mean(rhox0,'all').*g0.*ZZ;  Pl = Pt;  Pto = Pt; Ptoo = Pt; dPtdt = 0*Pt; dPtdto = dPtdt; dPtdtoo = dPtdto;
rhox   = rhox0.*(1+bPx.*(Pt-Pref));
rhom   = rhom0.*(1+bPm.*(Pt-Pref));
rho    = rhox;
M      = 0*rhox;
rhow   = (rho(icz(1:end-1),:)+rho(icz(2:end),:))/2;
rhou   = (rho(:,icx(1:end-1))+rho(:,icx(2:end)))/2;
rhoWo  = rhow.*W(:,2:end-1); rhoWoo = rhoWo; advn_mz = 0.*rhoWo(2:end-1,:);
rhoUo  = rhou.*U(2:end-1,:); rhoUoo = rhoUo; advn_mx = 0.*rhoUo;
mq = zeros(size(Tp));  xq = 1-mq;  
cmq    = c; cxq = c;  
cm_oxd = reshape(reshape(c,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
cm_oxd_all(:,:,cal.ioxd) = cm_oxd;
aT   = aTm;
kT   = kTm;
cP   = cPm; RhoCp = rho.*cP;
Adbt = aT./RhoCp;
Tp   = (Tp+273.15); %T = Tp;
T    = Tp.*exp(Adbt.*(Pt-Pref));
sm   = cPm.*log(Tp./Tref) + Dsm;  sx = cPx.*log(Tp./Tref);  
x    = xq;  m = mq; mu = m; chi = x;
dto  = dt;

% get volume fractions and bulk density
step    = 0;
EQtime  = 0;
FMtime  = 0;
TCtime  = 0;
UDtime  = 0;
a1      = 1; a2 = 0; a3 = 0; b1 = 1; b2 = 0; b3 = 0;

res  = 1;  tol = 1e-9;  it = 1; iter = 1;

while res > tol
    Pti = Pt; Ti = T; xi = xq;          

    % %%%%The melt model %%%%

        eqtime = tic;

        var.c      = reshape(c,Nx*Nz,cal.ncmp);   % component fractions [wt]
        var.T      = reshape(T,Nx*Nz,1)-273.15;   % temperature [C]
        var.P      = reshape(Pt,Nx*Nz,1)/1e9;     % pressure [GPa]
        var.m      = reshape(mq,Nx*Nz,1);         % melt fraction [wt](melt model 
        var.H2O    = var.c(:,end);                % water concentration [wt]
        var.X      = reshape(cm_oxd_all,Nz*Nx,9); % melt oxide fractions [wt %]
        cal.H2Osat = fluidsat(var); % water saturation [wt]

        %[var,cal] = meltmodel(var,cal,'E');
        [var,cal] = leappart(var,cal,'E');


        Tsol   = reshape(cal.Tsol,Nz,Nx);
        Tliq   = reshape(cal.Tliq,Nz,Nx);
        H2Osat = reshape(cal.H2Osat,Nz,Nx);

        mq = reshape(var.m.*(var.m>eps^0.5),Nz,Nx);
        xq = reshape(var.x.*(var.x>eps^0.5),Nz,Nx);
        mq = mq./(mq+xq);
        xq = xq./(mq+xq);
        x  = xq;  m = mq;

        cxq = reshape(var.cx,Nz,Nx,cal.ncmp);
        cmq = reshape(var.cm,Nz,Nx,cal.ncmp);
        cm  = cmq; cx = cxq;

        sref = 0;
        sm   = cPm.*log(Tp./Tref) + Dsm;  
        sx   = cPx.*log(Tp./Tref);

        eqtime = toc(eqtime);
        EQtime = EQtime + eqtime;

      if bndmode == 0 % Mid ocean Ridge set up 
        % Removing melt to get a suitable initial melt fraction
        if it>10 && any(m(:)>minit)
            m = m * (minit./max(m(:)))^0.25;
            SUM = x+m;
            x = x./SUM;  m = m./SUM;
            c = x.*cx + m.*cm;
            s = x.*sx + m.*sm;
            rho = 1./(m./rhom  + x./rhox);
        else
            s = x.*sx + m.*sm;
        end

      else
          s = x.*sx + m.*sm;
      end

        X    = rho.*x;
        M    = rho.*m;  RHO = X+M;
        C    = rho.*c;
        S    = rho.*s;

        update;
        Pf(2:end-1,2:end-1) = Pt;
        Px = Pt;

        [Tp,~ ] = StoT(Tp,S./rho,cat(3,Pt,Ptx)*0+Pref,cat(3,m,x),[cPm;cPx],[aTm;aTx],[bPm;bPx],cat(3,rhom0,rhox0),[sref+Dsm;sref],Tref,Pref);
        [T ,si] = StoT(T ,S./rho,cat(3,Pt,Ptx)       ,cat(3,m,x),[cPm;cPx],[aTm;aTx],[bPm;bPx],cat(3,rhom0,rhox0),[sref+Dsm;sref],Tref,Pref);
        sm = si(:,:,1); sx = si(:,:,2);

        res  = norm(Pt(:)-Pti(:),2)./norm(Pt(:),2) ...
             + norm( T(:)-Ti  (:),2)./norm( T(:),2);

        it = it+1;

end

Pto  = Pt;
To   = T;  
Tpo  = Tp;
So   = S;  
Mo   = M;
Co   = C;
Xo   = X;
rhoo = rho;

sm   = cPm.*log(Tp./Tref) + Dsm;  sx = cPx.*log(Tp./Tref);  


% get trace element phase compositions
Ktrc = zeros(Nz,Nx,cal.ntrc);
trcm = zeros(Nz,Nx,cal.ntrc);
trcx = zeros(Nz,Nx,cal.ntrc);
for i = 1:cal.ntrc
    for j=1:cal.nmem; Ktrc(:,:,i) = Ktrc(:,:,i) + cal.Ktrc_mem(i,j) .* c_mem(:,:,j)./100; end

    trcm(:,:,i)  = trc(:,:,i)./(m + x.*Ktrc(:,:,i));
    trcx(:,:,i)  = trc(:,:,i)./(m./Ktrc(:,:,i) + x);
end

% get geochemical component densities
TRC = zeros(Nz,Nx,cal.ntrc);
for i = 1:cal.ntrc
    TRC(:,:,i)  = rho.*(m.*trcm(:,:,i) + x.*trcx(:,:,i));
end
TRCo = TRC;

% initialise phase change rates
Gx  = 0.*x; Gm  = 0.*m; 
Gem = 0.*m; Gex = 0.*x; 
Gemc = 0.*c; Gexc = 0.*c;
Gemt = 0.*trc; Gext = 0.*trc;

% initialise auxiliary variables
dSdt   = 0.*T;  dSdto  = dSdt; diss_h = 0.*T;
dTpdt  = 0.*T;  dTpdto = dTpdt;
dTdt   = 0.*T;  dTdto  = dTdt;
dCdt   = 0.*c;  dCdto  = dCdt;
dXdt   = 0.*x;  dXdto  = dXdt;
dMdt   = 0.*m;  dMdto  = dMdt;
bnd_TRC = zeros(Nz,Nx,cal.ntrc);
adv_TRC = zeros(Nz,Nx,cal.ntrc);
dff_TRC = zeros(Nz,Nx,cal.ntrc);
K_trc     = zeros(Nz,Nx,cal.ntrc);
dTRCdt  = 0.*trc; dTRCdto = dTRCdt;
upd_S   = 0.*S;
upd_Tp  = 0.*Tp;
upd_T   = 0.*T;
upd_C   = 0.*C;
upd_X   = 0.*X;
upd_M   = 0.*M;
upd_rho = 0.*rho;
upd_eta = 0.*eta;
upd_TRC = 0.*TRC;
% upd_IR = 0.*IR;

% initialise timing and iterative parameters
frst    = 1;
step    = 0;
time    = 0;
iter    = 2;
hist    = [];
dsumSdto = 0;  dsumSdt = 0;
dsumBdto = 0;  dsumBdt = 0;
dsumMdto = 0;  dsumMdt = 0;
dsumXdto = 0;  dsumXdt = 0;
dsumCdto = 0;  dsumCdt = 0;
dsumTdto = 0;  dsumTdt = 0;

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = [outdir,'/',runID,'/',runID,'_cont.mat'];
    elseif restart > 0  % restart from specified continuation frame
        name = [outdir,'/',runID,'/',runID,'_',num2str(restart),'.mat'];
    end
    
    if exist(name,'file')
        fprintf('\n   restart from %s \n\n',name);
        load(name,'U','W','Pf','Pc','Pt','x','m','xq','mq','chi','mu','X','M','S','C','T','Tp','c','cm','cx','TRC','trc','dSdt','dCdt','dXdt','dMdt','drhodt','dTRCdt','Gx','Gm','Gem','Gex','rho','eta','eII','tII','dt','time','step','VolSrc','Div_V','qDz','qDx','wx','wm','cal');
        name = [outdir,'/',runID,'/',runID,'_hist'];
        load(name,'hist');

        SOL = [W(:);U(:);Pf(:);Pc(:)];
        RHO = X+M; 
     
        update;
        phseql;
        store;
        fluidmech;
        update;
        output;

        time    = time+dt;
        step    = step+1;

    else % continuation file does not exist, start from scratch
        fprintf('\n   !!! restart file does not exist !!! \n   => starting run from scratch %s \n\n',runID);
        store;
        fluidmech;
        update;
        history;
        output;
    end
else
    % complete, plot, and save initial condition
    store;
    fluidmech;
    update;
    history;
    output;
    step = step+1;
end

restart = 0;
