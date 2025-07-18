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
cx_msy = reshape(reshape(cx_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);

% update mineral systems oxide compositions for solid assemblage
cx_msy_oxd = zeros(Nz,Nx,cal.nmsy,cal.noxd);
for j = 1:cal.nmsy
    cx_msy_oxd(:,:,j,:) = reshape(reshape(cx_mem(:,:,cal.msy_mem(j,:)==1),Nz*Nx,sum(cal.msy_mem(j,:)==1))*cal.mem_oxd(cal.msy_mem(j,:)==1,:)./sum(reshape(cx_mem(:,:,cal.msy_mem(j,:)==1),Nz*Nx,sum(cal.msy_mem(j,:)==1))+1e-32,2),Nz,Nx,1,cal.noxd);
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

% interpolate to staggered stencil nodes
rhoxw  = (rhox(icz(1:end-1),:)+rhox(icz(2:end),:))/2;
rhomw  = (rhom(icz(1:end-1),:)+rhom(icz(2:end),:))/2;

rhow   = (rho(icz(1:end-1),:)+rho(icz(2:end),:))/2;
rhou   = (rho(:,icx(1:end-1))+rho(:,icx(2:end)))/2;

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

% interpolate to staggered stencil nodes
muw  = (mu (icz(1:end-1),icx)+mu (icz(2:end),icx))./2;
muu  = (mu (icz,icx(1:end-1))+mu (icz,icx(2:end)))./2;

chi_mem = reshape(reshape(cx_mem/100.*rhox,Nz*Nx,cal.nmem)./cal.rhox0,Nz,Nx,cal.nmem);
chi_mem = chi_mem./sum(chi_mem,3);

% update thermal parameters
aT = mu.*aTm + chi.*aTx;
bP = mu.*bPm + chi.*bPx;
kT = mu.*kTm + chi.*kTx;
cP = mu.*cPm + chi.*cPx;
RhoCp = mu.*rhom.*cPm + chi.*rhox.*cPx;
Adbt  = mu.*aTm./rhom./cPm + chi.*aTx./rhox./cPx;

% Extracted bounday conditions
if iter==1 % update two-phase masking once per time step
    twophs  = double(mu (icz,icx)>=mulim); 
    twophsw  = (twophs (icz(1:end-1),icx)+twophs (icz(2:end),icx))./2;
    twophsu  = (twophs (icz,icx(1:end-1))+twophs (icz,icx(2:end)))./2;
    % twophsw = double(muw         >=mulim);
    % twophsu = double(muu         >=mulim);
end  

% update lithostatic pressure
if Nz==1; Pt    = max(1e7,(1-alpha).*Pt + alpha.*(Ptop.*ones(size(Tp)) + Pcouple*(Pchmb + Pf(2:end-1,2:end-1)))); else
    Pl(1,:)     = repmat(mean(rhow(1,:),2).*g0.*h/2,1,Nx) + Ptop;
    Pl(2:end,:) = Pl(1,:) + repmat(cumsum(mean(rhow(2:end-1,:),2).*g0.*h),1,Nx);
    Pt          = max(1e7,(1-alpha).*Pt + alpha.*(Pl + Pcouple*(Pchmb + Pf(2:end-1,2:end-1))));
end
Ptx = Pt + Pcouple.*Pc(2:end-1,2:end-1)./(1-mu);

% update pure phase viscosities
etam   = reshape(Giordano08(reshape(cm_oxd_all,Nz*Nx,9),T(:)-273.15),Nz,Nx);  % T in [C]
etax0  = reshape(prod(cal.etax0(1:end-1).^reshape(chi_mem(:,:,1:end-1)+eps,Nz*Nx,cal.nmem-1),2),Nz,Nx);
etax   = etax0 .* exp(cal.Eax./(8.3145.*T)-cal.Eax./(8.3145.*(T1+273.15)));

% get coefficient contrasts
kv = permute(cat(3,etax,etam),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,2),[4,1,2,3])./permute(repmat(kv,1,1,1,2),[1,4,2,3]);

% get permission weights
ff = permute(cat(3,chi,mu+mulim),[3,1,2]);
FF = permute(repmat(ff,1,1,1,2),[4,1,2,3]);
Sf = (FF./cal.BB).^(1./cal.CC);  Sf = Sf./sum(Sf,2);
Xf = sum(cal.AA.*Sf,2).*FF + (1-sum(cal.AA.*Sf,2)).*Sf;

% get momentum flux and transfer coefficients
thtv = squeeze(prod(Mv.^Xf,2));
Kv   = ff.*kv.*thtv;
Cv   = Kv./dx0.^2;

% get effective viscosity
eta0  = squeeze(sum(Kv,1)); if Nx==1; eta0 = eta0.'; end

% get yield viscosity
etay   = tyield./(eII + eps^1.25) + etaymin;
eta    = eta.^gamma .* ((1./etay + 1./eta0).^-1).^(1-gamma);

% traditional two-phase coefficients
KD     = (mu+mulim).^2./squeeze(Cv(2,:,:));  % melt segregation coeff
KD     = (1./KD + 1./1e-10).^-1 + 1e-24;

zeta0  = eta./(mu+mulim);  % solid compaction coeff

% get yield viscosity
zetay  = (1-twophs(2:end-1,2:end-1)).*pyield/eps^1.25 + twophs(2:end-1,2:end-1).*pyield./(max(0,Div_V)+eps^1.25) + etaymin./(mu+mulim);
zeta   = zeta.*gamma + ((1./zetay + 1./zeta0).^-1).*(1-gamma);

% interpolate to staggered stencil nodes
etaco  = (eta(icz(1:end-1),icx(1:end-1)).*eta(icz(2:end),icx(1:end-1)) ...
       .* eta(icz(1:end-1),icx(2:end  )).*eta(icz(2:end),icx(2:end  ))).^0.25;

KDw    = (KD(icz(1:end-1),:) + KD(icz(2:end),:)).*0.5;
KDu    = (KD(:,icx(1:end-1)) + KD(:,icx(2:end))).*0.5;
% KDw    = (KD(icz(1:end-1),:) .* KD(icz(2:end),:)).^0.5;
% KDu    = (KD(:,icx(1:end-1)) .* KD(:,icx(2:end))).^0.5;

% update velocity magnitudes
Vel = sqrt(((W(1:end-1,2:end-1)+W(2:end,2:end-1))/2).^2 ...
         + ((U(2:end-1,1:end-1)+U(2:end-1,2:end))/2).^2);

qD  = sqrt(((qDz(1:end-1,2:end-1)+qDz(2:end,2:end-1))/2).^2 ...
         + ((qDx(2:end-1,1:end-1)+qDx(2:end-1,2:end))/2).^2);

% update velocity divergences
Div_V  = ddz(W  (:,2:end-1),h) + ddx(U  (2:end-1,:),h);                    % get velocity divergence
Div_qD = ddz(qDz(:,2:end-1),h) + ddx(qDx(2:end-1,:),h);                    % get velocity divergence

% update strain rates
exx = diff(U(2:end-1,:),1,2)./h - Div_V./3;                                % x-normal strain rate
ezz = diff(W(:,2:end-1),1,1)./h - Div_V./3;                                % z-normal strain rate
exz = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h);                                % shear strain rate

eII = (0.5.*(exx.^2 + ezz.^2 ...
       + 2.*(exz(1:end-1,1:end-1).^2+exz(2:end,1:end-1).^2 ...
       +     exz(1:end-1,2:end  ).^2+exz(2:end,2:end  ).^2)/4)).^0.5 + eps;% strain rate magnitude

etamax = etacntr.*max(min(eta(:)),etamin);
eta    = 1./(1./etamax + 1./eta) + etamin;
zetamax = 1./(1./(eta./mulim) + 1./(etamax./(mu+mulim)));
zetamin = etamin./(mu+mulim);
zeta   = 1./(1./zetamax + 1./zeta) + zetamin;

% update dimensionless numbers
Ra     = Vel.*D/10./((kT)./rho./cP);
Re     = Vel.*D/10./( eta./rho    );
Rum    = abs(wm(1:end-1,2:end-1)+wm(2:end,2:end-1))/2./Vel;
delta  = sqrt(KD.*zeta);

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
    diss = kT./T.*(grdTz (2:end-1,2:end-1).^2 + grdTx (2:end-1,2:end-1).^2) ...
         + exx.*txx + ezz.*tzz ...
         + 2.*(exz(1:end-1,1:end-1)+exz(2:end,1:end-1)+exz(1:end-1,2:end)+exz(2:end,2:end))./4 ...
            .*(txz(1:end-1,1:end-1)+txz(2:end,1:end-1)+txz(1:end-1,2:end)+txz(2:end,2:end))./4 ...
         +  KD .* ((qDz(1:end-1,2:end-1)+qDz(2:end,2:end-1))./2).^2 ...
         +  KD .* ((qDx(2:end-1,1:end-1)+qDx(2:end-1,2:end))./2).^2;
end

% update time step
dtk = (h/2)^2/max(kT(:)./rho(:)./cP(:));                       % diffusive time step size
dta =  h/2   /max(abs([Um(:);Wm(:);Ux(:);Wx(:)]));             % advective time step size
dt  = min([1.1*dto,min(CFL*[dtk,dta]),dtmax]);                 % time step size

UDtime = UDtime + toc;
