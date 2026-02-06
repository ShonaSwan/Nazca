%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************
tic;

% update phase oxide compositions

c_oxd  = reshape(reshape(c ,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
cm_oxd = reshape(reshape(cm,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
cx_oxd = reshape(reshape(cx,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);

% update phase mineral end-member compositions
c_mem  = reshape(reshape(c ,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);
cm_mem = reshape(reshape(cm,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);
cx_mem = reshape(reshape(cx,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);

% update mineral systems composition for solid assemblage
c_msy  = reshape(reshape( c_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);
cm_msy = reshape(reshape(cm_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);
cx_msy = reshape(reshape(cx_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);

% update mineral systems oxide compositions for solid assemblage
cx_msy_oxd = zeros(Nz,Nx,cal.nmsy,cal.noxd);
for j = 1:cal.nmsy
    cx_msy_oxd(:,:,j,:) = reshape(reshape(cx_mem(:,:,cal.msy_mem(j,:)==1),Nz*Nx,nnz(cal.msy_mem(j,:)==1))*cal.mem_oxd(cal.msy_mem(j,:)==1,:)./sum(reshape(cx_mem(:,:,cal.msy_mem(j,:)==1),Nz*Nx,nnz(cal.msy_mem(j,:)==1))+1e-32,2),Nz,Nx,1,cal.noxd);
end

cm_oxd_all = zeros(size(c,1),size(c,2),9);
cm_oxd_all(:,:,cal.ioxd) = cm_oxd;
cx_oxd_all = zeros(size(c,1),size(c,2),9);
cx_oxd_all(:,:,cal.ioxd) = cx_oxd;
 c_oxd_all = zeros(size(c,1),size(c,2),9);
 c_oxd_all(:,:,cal.ioxd) = c_oxd;

% update phase densities
rhom0  = reshape(DensityX(reshape(cm_oxd_all,Nz*Nx,9),293,Pref./1e8)    ,Nz,Nx);
rhox0  = reshape(sum(reshape(cx_mem/100,Nz*Nx,cal.nmem)./cal.rhox0,2).^-1,Nz,Nx);

rhom   = rhom0 .* (1 - aTm.*(T-293) + bPm.*(Pt-Pref));
rhox   = rhox0 .* (1 - aTx.*(T-293) + bPx.*(Pt-Pref));

rho0   = 1./(m./rhom0 + x./rhox0);
rho    = 1./(m./rhom  + x./rhox );

% S   = rho.*s;
% C   = rho.*c;
% M   = rho.*m;
% X   = rho.*x;
% TRC = rho.*trc;

% interpolate to staggered stencil nodes

if meansw == 0 % Geometric
   rhoxw  = (rhox(icz(1:end-1),:).*rhox(icz(2:end),:)).^0.5; 
   rhoxw(end,:) = rhoxw(end-1,:).^2 ./ rhoxw(end-2,:);
   rhoxu  = (rhox(:,icx(1:end-1)).*rhox(:,icx(2:end))).^0.5;
   rhomw  = (rhom(icz(1:end-1),:).*rhom(icz(2:end),:)).^0.5; 
   rhomu  = (rhom(:,icx(1:end-1)).*rhom(:,icx(2:end))).^0.5;
elseif meansw == 1 % Arithmetic
   rhoxw  = (rhox(icz(1:end-1),:)+rhox(icz(2:end),:))/2;   
   rhoxw(end,:) = 2*rhoxw(end-1,:) - rhoxw(end-2,:);
   rhoxu  = (rhox(:,icx(1:end-1))+rhox(:,icx(2:end)))/2;
   rhomw  = (rhom(icz(1:end-1),:)+rhom(icz(2:end),:))/2; 
   rhomu  = (rhom(:,icx(1:end-1))+rhom(:,icx(2:end)))/2;  
end


if meansw == 0 % Geometric
   rhow   = (rho(icz(1:end-1),:).*rho(icz(2:end),:)).^0.5; 
   rhow(end,:) = rhow(end-1,:).^2 ./ rhow(end-2,:);
   rhou   = (rho(:,icx(1:end-1)).*rho(:,icx(2:end))).^0.5;   
elseif meansw == 1 % Arithmetic
   rhow   = (rho(icz(1:end-1),:)+rho(icz(2:end),:))/2;   
   rhow(end,:) = 2*rhow(end-1,:) - rhow(end-2,:);
   rhou   = (rho(:,icx(1:end-1))+rho(:,icx(2:end)))/2;          
end


if meansw == 0 % Geometric
   Mw     = (M(icz(1:end-1),:).*M(icz(2:end),:)).^0.5;       
   Mu     = (M(:,icx(1:end-1)).*M(:,icx(2:end))).^0.5; 
   Xw     = (X(icz(1:end-1),:).*X(icz(2:end),:)).^0.5;
   Xu     = (X(:,icx(1:end-1)).*X(:,icx(2:end))).^0.5;
elseif meansw == 1 % Arithmetic
   Mw     = (M(icz(1:end-1),:)+M(icz(2:end),:))/2;             
   Mu     = (M(:,icx(1:end-1))+M(:,icx(2:end)))/2;   
   Xw     = (X(icz(1:end-1),:)+X(icz(2:end),:))/2;
   Xu     = (X(:,icx(1:end-1))+X(:,icx(2:end)))/2;
end
 

if meansw == 0 % Geometric
   mz     = (m(icz(1:end-1),:).*m(icz(2:end),:)).^0.5;        
   mx     = (m(:,icx(1:end-1)).*m(:,icx(2:end))).^0.5;       
elseif meansw == 1 % Arithmetic
   mz     = (m(icz(1:end-1),:)+m(icz(2:end),:))/2;      
   mx     = (m(:,icx(1:end-1))+m(:,icx(2:end)))/2;                  
end


% update density contrasts
Drhow  = rhow -mean(rhow,2);
Drhomw = rhomw-mean(rhow,2);

% extract potential density
rhomp = rhom0 .* (1 - aTm.*(Tp-293));
rhoxp = rhox0 .* (1 - aTx.*(Tp-293));

rhop  = 1./(m./rhomp + x./rhoxp);

% convert weight to volume fraction, update bulk density
chi    = max(0,min(1, x.*rho./rhox ));
mu     = max(0,min(1, m.*rho./rhom ));

mucff  = (1./mu.^4 + 1./mumax.^4).^-(1/4) + mumin;

% interpolate to staggered stencil nodes
 
if meansw == 0 % Geometric
   muw  = (mu (icz(1:end-1),:).*mu (icz(2:end),:)).^0.5; 
   muu  = (mu (:,icx(1:end-1)).*mu (:,icx(2:end))).^0.5;     
   chiw = (chi(icz(1:end-1),:).*chi(icz(2:end),:)).^0.5;
   chiu = (chi(:,icx(1:end-1)).*chi(:,icx(2:end))).^0.5;    
elseif meansw == 1 % Arithmetic
   muw  = (mu (icz(1:end-1),:)+mu (icz(2:end),:))./2;    
   muu  = (mu (:,icx(1:end-1))+mu (:,icx(2:end)))./2;  
   chiw = (chi (icz(1:end-1),:)+chi (icz(2:end),:))./2;
   chiu = (chi (:,icx(1:end-1))+chi (:,icx(2:end)))./2;
end

chi_mem = reshape(reshape(cx_mem/100.*rhox,Nz*Nx,cal.nmem)./cal.rhox0,Nz,Nx,cal.nmem);
chi_mem = chi_mem./sum(chi_mem,3);

% update thermal parameters
aT    = mu.*aTm + chi.*aTx;
bP    = mu.*bPm + chi.*bPx;
kT    = mu.*kTm + chi.*kTx;
cP    = mu.*cPm + chi.*cPx;
RhoCp = mu.*rhom.*cPm + chi.*rhox.*cPx;
Adbt  = mu.*aTm./rhom./cPm + chi.*aTx./rhox./cPx;

% Extracted bounday conditions
if iter<3 % update two-phase masking once per time step 
    twophs  = double(mu (icz,icx)>mumin);
    twophsw = double(muw(:  ,icx)>mumin);
    twophsu = double(muu(icz,:  )>mumin);
end 

% update lithostatic pressure
if Nz==1; Pt    = max(Ptop,(1-delta).*Pt + delta.*(Ptop.*ones(size(Tp)) + Pcouple*Pf(2:end-1,2:end-1))); else
    Pl(1,:)     = repmat(mean(rhow(1,:),2).*g0.*h/2,1,Nx) + Ptop;
    Pl(2:end,:) = Pl(1,:) + repmat(cumsum(mean(rhow(2:end-1,:),2).*g0.*h),1,Nx);
    Pt          = max(Ptop,(1-delta).*Pt + delta.*(Pl + Pcouple*Pf(2:end-1,2:end-1)));   
end
Ptx = Pt + Pcouple.*Pc(2:end-1,2:end-1)./(1-mucff);

% update melt viscosity
etam   = reshape(Giordano08(reshape(cm_oxd_all,Nz*Nx,9),T(:)-273.15),Nz,Nx);  % T in [C]
etamax = etacntr.*min(etam(:));
etam   = 1./(1./etamax + 1./etam);

% update solid T,C-dependent viscosity
etax0  = reshape(prod(cal.etax0(1:end-1).^reshape(chi_mem(:,:,1:end-1)+eps,Nz*Nx,cal.nmem-1),2),Nz,Nx);
etax0  = etax0 .* exp(cal.Eax./(8.3145.*T)-cal.Eax./(8.3145.*(T1+273.15)));

% diffusion creep viscosity (grain size dependence)
dx_ref   = 0.001;  
eta_diff = etax0 .* (dx0./dx_ref).^2;

% dislocation creep viscosity (strain rate dependence)
eps_ref  = 1e-14;  
eta_disl = etax0 .* (max(eII, 1e-18) / eps_ref).^((1/n_disl)-1);

% composite matrix rheology
etax     = (1./eta_diff + 1./eta_disl).^-1;

% effective two-phase coefficients
eta0   = etax.*exp(-lmbd_melt.*mucff); % matrix shear viscosity
zeta0  = eta0./mucff;                  % matrix compaction viscosity
kphi   = dx0^2./b_perm .* mucff.^3;    % matrix permeability
KD     = kphi./etam;% + (h/100)^2./zeta; % Darcy coefficient
% KDi    = (1-delta).*KDi + delta.*KD0;%./(1 + max(0,(wm(1:end-1,2:end-1)+wm(2:end,2:end-1))/2)./rms(wm(:)+eps));
Ks     = KD./mucff;                    % segregation drag coefficient

% apply min/max bounds to viscosities
etamax  = etacntr.*max(min(eta0(:)),etamin);
eta0    = 1./(1./etamax + 1./eta0) + etamin;
zetamax = etamax./mucff;
zetamin = etamin./mucff;
zeta0   = 1./(1./zetamax + 1./zeta0) + zetamin;

% get shear and compaction strain rates

% update velocity divergences
Div_V    = ddz(W   (:,2:end-1),h) + ddx(U   (2:end-1,:),h); % get matrix velocity divergence
Div_DV   = ddz(wm  (:,2:end-1),h) + ddx(um  (2:end-1,:),h); % get segregation velocity divergence
Div_Vmix = ddz(Wmix(:,2:end-1),h) + ddx(Umix(2:end-1,:),h); % get mixture velocity divergence
rhoW     = rhow(:,icx).*W;
rhoU     = rhou(icz,:).*U;
Mwm      = Mw(:,icx).*wm;
Mum      = Mu(icz,:).*um;
Xwx      = Xw(:,icx).*W;
Xux      = Xu(icz,:).*U;
Div_rhoV = ddz(rhoW(:,2:end-1),h) + ddx(rhoU(2:end-1,:),h) ...
         + ddz(Mwm (:,2:end-1),h) + ddx(Mum (2:end-1,:),h);
DivRhoxV = ddz(rhoxw.*W(:,2:end-1),h) + ddx(rhoxu.*U(2:end-1,:),h);

% update compaction rate
if step>0 && ~restart
    % ups  = - 1./chi.*((a1*chi-a2*chio-a3*chioo)/dt - advect(chi,U(2:end-1,:),W(:,2:end-1),h,{ADVN,'vdf'},[1,2],BCA) - Gx./rhox);
    ups = ((a1*rhox-a2*rhoxo-a3*rhoxoo)/dt + DivRhoxV - Gex)./rhox;
end

% update strain rates
exx = diff(U(2:end-1,:),1,2)./h - Div_V./3;                                % x-normal strain rate
ezz = diff(W(:,2:end-1),1,1)./h - Div_V./3;                                % z-normal strain rate
exz = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h);                                % shear strain rate

eII = (0.5.*(exx.^2 + ezz.^2 ...
       + 2.*(exz(1:end-1,1:end-1).^2+exz(2:end,1:end-1).^2 ...
       +     exz(1:end-1,2:end  ).^2+exz(2:end,2:end  ).^2)/4)).^0.5 + eps;

mufact = (1-exp(-mu./mumin));
Peff   = Pt - mufact.*(Pl + Pf(2:end-1,2:end-1));
yieldt = max(1e4,mufact.*min(pyield + Pc(2:end-1,2:end-1), 100*pyield - Pc(2:end-1,2:end-1)) + (1-mufact).*(tyield + 0.5*Pt));
yieldp = max(1e4,pyield - tII);

% get yield shear viscosity
etai   = (1-delta).*etai + delta.*(yieldt./(eII + 1e-32) + etaymin);
eta    = ((1./etai.^2 + 1./eta0.^2).^-(1/2));
zeta0  = zeta0.*min(1,etai./eta0);

% get yield compaction viscosity
upsy   = mufact.*(max(0,ups)+max(0,-ups/100));
zetai  = (1-delta).*zetai + delta.*(yieldp./(upsy + 1e-32) + etaymin);
zeta   = ((1./zetai.^2 + 1./zeta0.^2).^-(1/2));

if cff_reg
    eta = log10(eta);
    for i = 1:cff_reg
        eta = eta + diffus(eta,1/8*ones(size(eta)),1,[1,2],BCD);
    end
    eta = 10.^eta;

    zeta = log10(zeta);
    for i = 1:cff_reg
        zeta = zeta + diffus(zeta,1/8*ones(size(zeta)),1,[1,2],BCD);
    end
    zeta = 10.^zeta;

    KD = log10(KD);
    for i = 1:cff_reg
        KD = KD + diffus(KD,1/8*ones(size(KD)),1,[1,2],BCD);
    end
    KD = 10.^KD;
    
    Ks = log10(Ks);
    for i = 1:cff_reg
        Ks = Ks + diffus(Ks,1/8*ones(size(Ks)),1,[1,2],BCD);
    end
    Ks = 10.^Ks;
% else
%      eta =  etai;
%     zeta = zetai;
%     KD   = KDi;
%     Ks   = Ksi;
end

% interpolate to staggered stencil nodes
if     meansw == 0 % Geometric
    etaco  = (eta(icz(1:end-1),icx(1:end-1)).*eta(icz(2:end),icx(1:end-1)) ...
           .* eta(icz(1:end-1),icx(2:end  )).*eta(icz(2:end),icx(2:end  ))).^0.25;
elseif meansw == 1 % Arithmetic
    etaco  = (eta(icz(1:end-1),icx(1:end-1))+eta(icz(2:end),icx(1:end-1)) ...
           +  eta(icz(1:end-1),icx(2:end  ))+eta(icz(2:end),icx(2:end  ))).*0.25;
end

if meansw == 0 % Geometric
   Ksw    = (Ks(icz(1:end-1),:) .* Ks(icz(2:end),:)).^0.5;  
   Ksu    = (Ks(:,icx(1:end-1)) .* Ks(:,icx(2:end))).^0.5;        
elseif meansw == 1 % Arithmetic
   Ksw    = (Ks(icz(1:end-1),:) + Ks(icz(2:end),:)).*0.5; 
   Ksu    = (Ks(:,icx(1:end-1)) + Ks(:,icx(2:end))).*0.5;                   
end

% update velocity magnitudes
Vx  = sqrt(((W(1:end-1,2:end-1)+W(2:end,2:end-1))/2).^2 ...
         + ((U(2:end-1,1:end-1)+U(2:end-1,2:end))/2).^2);
Vm  = sqrt(((wm(1:end-1,2:end-1)+wm(2:end,2:end-1))/2).^2 ...
         + ((um(2:end-1,1:end-1)+um(2:end-1,2:end))/2).^2);
qD  = sqrt(((qDz(1:end-1,2:end-1)+qDz(2:end,2:end-1))/2).^2 ...
         + ((qDx(2:end-1,1:end-1)+qDx(2:end-1,2:end))/2).^2);

% update melt dispersivity
kd = Vm.*Delta + km;

% update dimensionless numbers
Ra     = Vx.*D/10./((kT)./rho./cP);
Re     = Vx.*D/10./( eta./rho    );
Rc     = Vx./Vm;
deltac = sqrt(KD.*zeta);

% update stresses
txx = eta   .* exx;                                                        % x-normal stress
tzz = eta   .* ezz;                                                        % z-normal stress
txz = etaco .* exz;                                                        % xz-shear stress

tII = (0.5.*(txx.^2 + tzz.^2 ...
       + 2.*(txz(1:end-1,1:end-1).^2+txz(2:end,1:end-1).^2 ...
       +     txz(1:end-1,2:end  ).^2+txz(2:end,2:end  ).^2)/4)).^0.5 + eps;

% heat dissipation (entropy production) rate
if Nz==1 && Nx==1
    diss = 0.*T;  % no dissipation in 0-D mode (no diffusion, no shear deformation, no segregation)
else
    [grdTx ,grdTz ] = gradient(T(icz,icx),h);
    exz_ce = (exz(1:end-1,1:end-1)+exz(2:end,1:end-1)+exz(1:end-1,2:end)+exz(2:end,2:end))./4;
    diss = min(1e-3,diss ...
         + kT./T.*(grdTz (2:end-1,2:end-1).^2 + grdTx (2:end-1,2:end-1).^2) ...
         + eta.*exx.^2 + eta.*ezz.^2 + 2.*eta.*exz_ce.^2 ...
         + mucff./Ks .* ((wm(1:end-1,2:end-1)+wm(2:end,2:end-1))./2).^2 ...
         + mucff./Ks .* ((um(2:end-1,1:end-1)+um(2:end-1,2:end))./2).^2)/2;
end

% update time step
dtk = (h/2)^2/max([kT(:)./rho(:)./cP(:);kd(:)]);                           % diffusive time step size
dta =  h/2   /max(abs([Um(:);Wm(:);Ux(:);Wx(:)]));                         % advective time step size
dt  = min([2*dto,min(CFL*[dtk,dta]),dtmax]);     % time step size

UDtime = UDtime + toc;
