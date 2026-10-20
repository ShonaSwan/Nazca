%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************
tic;

% update phase oxide compositions

maj_oxd  = reshape(reshape(maj ,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
majm_oxd = reshape(reshape(majm,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
majx_oxd = reshape(reshape(majx,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);

% update phase mineral end-member compositions
maj_mem  = reshape(reshape(maj ,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);
majm_mem = reshape(reshape(majm,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);
majx_mem = reshape(reshape(majx,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);

% update mineral systems composition for solid assemblage
maj_msy  = reshape(reshape( maj_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);
majm_msy = reshape(reshape(majm_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);
majx_msy = reshape(reshape(majx_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);

% update mineral systems oxide compositions for solid assemblage
majx_msy_oxd = zeros(Nz,Nx,cal.nmsy,cal.noxd);
for j = 1:cal.nmsy
    majx_msy_oxd(:,:,j,:) = reshape(reshape(majx_mem(:,:,cal.msy_mem(j,:)==1),Nz*Nx,nnz(cal.msy_mem(j,:)==1))*cal.mem_oxd(cal.msy_mem(j,:)==1,:)./sum(reshape(majx_mem(:,:,cal.msy_mem(j,:)==1),Nz*Nx,nnz(cal.msy_mem(j,:)==1))+1e-32,2),Nz,Nx,1,cal.noxd);
end

majm_oxd_all = zeros(size(maj,1),size(maj,2),9);
majx_oxd_all = zeros(size(maj,1),size(maj,2),9);
 maj_oxd_all = zeros(size(maj,1),size(maj,2),9);
if cal.noxd>9
    majm_oxd_all = majm_oxd(:,:,cal.ioxd);
    majx_oxd_all = majx_oxd(:,:,cal.ioxd);
     maj_oxd_all =  maj_oxd(:,:,cal.ioxd);
else
    majm_oxd_all(:,:,cal.ioxd) = majm_oxd;
    majx_oxd_all(:,:,cal.ioxd) = majx_oxd;
     maj_oxd_all(:,:,cal.ioxd) = maj_oxd;
end

DivRhoxV = ddz(rhoxw.*W(:,2:end-1),h) + ddx(rhoxu.*U(2:end-1,:),h);

% update compaction rate
if step>0 && ~restart
    ups = twophs(2:end-1,2:end-1).*((a1*rhox-a2*rhoxo-a3*rhoxoo)/dt + DivRhoxV - Gex./chi)./rhox;
end

% update phase densities
rhom0  = reshape(DensityX(reshape(majm_oxd_all,Nz*Nx,9),293,Pref./1e8)    ,Nz,Nx);
rhox0  = reshape(sum(reshape(majx_mem/100,Nz*Nx,cal.nmem)./cal.rhox0,2).^-1,Nz,Nx);

rhom   = rhom0 .* (1 - aTm.*(T-293) + bPm.*(Pt-Pref));
rhox   = rhox0 .* (1 - aTx.*(T-293) + bPx.*(Pt-Pref));

rho0   = 1./(m./rhom0 + x./rhox0);
rho    = 1./(m./rhom  + x./rhox );

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

chi_mem = reshape(reshape(majx_mem/100.*rhox,Nz*Nx,cal.nmem)./cal.rhox0,Nz,Nx,cal.nmem);
chi_mem = chi_mem./sum(chi_mem,3);

% update thermal parameters
aT    = mu.*aTm + chi.*aTx;
bP    = mu.*bPm + chi.*bPx;
kT    = mu.*kTm + chi.*kTx;
cP    = mu.*cPm + chi.*cPx;
RhoCp = mu.*rhom.*cPm + chi.*rhox.*cPx;
Adbt  = mu.*aTm./rhom./cPm + chi.*aTx./rhox./cPx;

% update lithostatic pressure
if Nz==1; Pt    = max(Ptop,(1-delta).*Pt + delta.*(Ptop.*ones(size(Tp)) + Pcouple*Pf(2:end-1,2:end-1))); else
    Pl(1,:)     = repmat(mean(rhow(1,:),2).*g0.*h/2,1,Nx) + Ptop;
    Pl(2:end,:) = Pl(1,:) + repmat(cumsum(mean(rhow(2:end-1,:),2).*g0.*h),1,Nx);
    Pt          = max(Ptop,(1-delta).*Pt + delta.*(Pl + Pcouple*Pf(2:end-1,2:end-1)));   
end

% update melt viscosity
etam    = reshape(Giordano08(reshape(majm_oxd_all,Nz*Nx,9),T(:)-273.15),Nz,Nx);  % T in [C]
etammax = etacntr.*min(etam(:));
etam    = 1./(1./etammax + 1./etam);

% update solid T,C-dependent viscosity
etax0  = reshape(prod(cal.etax0.^reshape(chi_mem+eps,Nz*Nx,cal.nmem),2),Nz,Nx);
etax0  = etax0 .* exp(cal.Eax./(8.3145.*T)-cal.Eax./(8.3145.*(T1+273.15)));

% diffusion creep viscosity (grain size dependence)
dx_ref   = 0.001;  
eta_diff = etax0 .* (dx0./dx_ref).^2;

% dislocation creep viscosity (strain rate dependence)
eps_ref  = 1e-14;  
eta_disl = etax0 .* (max(eII, 1e-18) / eps_ref).^((1/n_disl)-1);

% composite matrix rheology
etax     = (1./eta_diff + 1./eta_disl).^-1;
eta0     = etax.*exp(-lmbd_melt.*mu) + etamin; % matrix shear viscosity

% get elastic moduli
G0  = reshape(prod(cal.G0.^reshape(chi_mem+eps,Nz*Nx,cal.nmem),2),Nz,Nx);
Gdt = G0.*dt;
Kdt = G0.*dt./sqrt(mu+1e-6);

if meansw == 0
    Gdtco = (Gdt(icz(1:end-1),icx(1:end-1)).*Gdt(icz(2:end),icx(1:end-1)) ...
          .* Gdt(icz(1:end-1),icx(2:end  )).*Gdt(icz(2:end),icx(2:end  ))).^0.25;
else
    Gdtco = (Gdt(icz(1:end-1),icx(1:end-1))+Gdt(icz(2:end),icx(1:end-1)) ...
          +  Gdt(icz(1:end-1),icx(2:end  ))+Gdt(icz(2:end),icx(2:end  ))).*0.25;
end

% effective two-phase coefficients
zeta0  = eta0./(mu+1e-6);              % matrix compaction viscosity
kphi   = dx0^2./b_perm .* (mu+1e-9).^3 .* (1-mu).^2;    % matrix permeability
KD     = kphi./etam;                   % Darcy coefficient
Ks     = KD./(mu+1e-9);                % segregation drag coefficient

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

% update strain rate components
exx = diff(U(2:end-1,:),1,2)./h - Div_V./3;                                % x-normal strain rate
ezz = diff(W(:,2:end-1),1,1)./h - Div_V./3;                                % z-normal strain rate
exz = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h);                                % shear strain rate

% update vorticity components
wxz = 1/2.*(diff(U,1,1)./h-diff(W,1,2)./h);
wzx = -wxz;

eII = (0.5.*(exx.^2 + ezz.^2 ...
       + 2.*(exz(1:end-1,1:end-1).^2+exz(2:end,1:end-1).^2 ...
       +     exz(1:end-1,2:end  ).^2+exz(2:end,2:end  ).^2)/4)).^0.5 + eps^2;

% update stresses
txx = eta   .* exx + tupd  .*(a2.*txxo + a3.*txxoo - advn_txx*dt);                                                        % x-normal stress
tzz = eta   .* ezz + tupd  .*(a2.*tzzo + a3.*tzzoo - advn_tzz*dt);                                                        % z-normal stress
txz = etaco .* exz + tupdco.*(a2.*txzo + a3.*txzoo - advn_txz*dt);                                                        % xz-shear stress

tII = (0.5.*(txx.^2 + tzz.^2 ...
       + 2.*(txz(1:end-1,1:end-1).^2+txz(2:end,1:end-1).^2 ...
       +     txz(1:end-1,2:end  ).^2+txz(2:end,2:end  ).^2)/4)).^0.5 + eps^2;

% get elastic strain rate components
exxe   = (a1*txx-a2*txxo-a3*txxoo+advn_txx*dt)./Gdt;
ezze   = (a1*tzz-a2*tzzo-a3*tzzoo+advn_tzz*dt)./Gdt;
exze   = (a1*txz-a2*txzo-a3*txzoo+advn_txz*dt)./Gdtco;

eIIe   = (0.5.*(exxe.^2 + ezze.^2 ...
       + 2.*(exze(1:end-1,1:end-1).^2+exze(2:end,1:end-1).^2 ...
       +     exze(1:end-1,2:end  ).^2+exze(2:end,2:end  ).^2)/4)).^0.5 + eps^2;

% get visco-plastic strain rate components
exxvp = txx./etavp;
ezzvp = tzz./etavp;
exzvp = txz./etavpco;

eIIvp  = (0.5.*(exxvp.^2 + ezzvp.^2 ...
       + 2.*(exzvp(1:end-1,1:end-1).^2+exzvp(2:end,1:end-1).^2 ...
       +     exzvp(1:end-1,2:end  ).^2+exzvp(2:end,2:end  ).^2)/4)).^0.5 + eps^2;

% get elastic and viscoplastic decompaction rate components
upse   = -(a1*Pc(2:end-1,2:end-1)-a2*Pco-a3*Pcoo+advn_Pc*dt)./Kdt;
upsvp  = -Pc(2:end-1,2:end-1)./zetavp;

% get stress advection with rotation using Jaumann derivative
wxzce  = (wxz(1:end-1,1:end-1)+wxz(2:end,1:end-1)+wxz(1:end-1,2:end)+wxz(2:end,2:end))/4;
txzce  = (txz(1:end-1,1:end-1)+txz(2:end,1:end-1)+txz(1:end-1,2:end)+txz(2:end,2:end))/4;
txxco  = (txx(icz(1:end-1),icx(1:end-1))+txx(icz(2:end),icx(1:end-1))+txx(icz(1:end-1),icx(2:end))+txx(icz(2:end),icx(2:end)))/4;
tzzco  = (tzz(icz(1:end-1),icx(1:end-1))+tzz(icz(2:end),icx(1:end-1))+tzz(icz(1:end-1),icx(2:end))+tzz(icz(2:end),icx(2:end)))/4;

advn_txx = advect(txx,U(2:end-1,:),W(:,2:end-1),h,{ADVN,'vdf'},[1,2],BCA) ...
         - 2.*txzce.*wxzce;
advn_tzz = advect(tzz,U(2:end-1,:),W(:,2:end-1),h,{ADVN,'vdf'},[1,2],BCA) ...
         + 2.*txzce.*wxzce;
advn_txz = advect(txz,(U(1:end-1,[1 1:end])+U(2:end,[1 1:end])+U(1:end-1,[1:end end])+U(2:end,[1:end end]))/4, ...
                      (W([1 1:end],1:end-1)+W([1 1:end],2:end)+W([1:end end],1:end-1)+W([1:end end],2:end))/4,h,{ADVN,'vdf'},[1,2],BCA) ...
         + (txxco-tzzco).*wxz;
advn_Pc  = advect(Pc(2:end-1,2:end-1),U(2:end-1,:),W(:,2:end-1),h,{ADVN,'vdf'},[1,2],BCA);

% get yield shear viscosities
if step>0
    mufact = twophs(2:end-1,2:end-1);%(1-exp(-mu./1e-4));
else
    mufact = 0;
end
frct   = reshape(sum(cal.frct0.*reshape(chi_mem+eps,Nz*Nx,cal.nmem),2),Nz,Nx);
if bndmode == 0
    frct = frct .* (1-0.5*exp(- XX     .^2./2/h^2 - ZZ.^2./2/h^2));
else
    frct = frct .* (1-0.5*exp(-(XX-L/2).^2./2/h^2 - ZZ.^2./2/h^2));
end
Peff   = Pl+P(2:end-1,2:end-1) - mufact.*(Pl + Pf(2:end-1,2:end-1));
yieldt = max(1e4,mufact.*min(pyield + Pc(2:end-1,2:end-1), 10*pyield - Pc(2:end-1,2:end-1)) + (1-mufact).*(cosd(frct).*tyield + sind(frct).*(Pl+P(2:end-1,2:end-1))));
yieldp = max(1e4,pyield - tII);

etay   = yieldt./eIIvp + etaymin;
% etay   = etay + (diff(etay(icz,:),2,1) + diff(etay(:,icx),2,2))/8;
etavp  = etavp.^(1-delta) .*  min(eta0,etay).^delta;

% get yield compaction viscosity
upsy   = mufact.*(max(eps^2,upsvp)+max(eps^2,-upsvp/10));
zetay  = yieldp./upsy + etaymin;
% zetay  = zetay + (diff(zetay(icz,:),2,1) + diff(zetay(:,icx),2,2))/8;
zetavp = zetavp.^(1-delta) .*  min(zeta0,zetay).^delta;

tupd = 1./(a1 + Gdt./ etavp);
pupd = 1./(a1 + Kdt./zetavp);

eta  = 1./(1./etavp  + a1./Gdt);
zeta = 1./(1./zetavp + a1./Kdt);


% interpolate to staggered stencil nodes
if     meansw == 0 % Geometric
    etavpco  = (etavp(icz(1:end-1),icx(1:end-1)).*etavp(icz(2:end),icx(1:end-1)) ...
             .* etavp(icz(1:end-1),icx(2:end  )).*etavp(icz(2:end),icx(2:end  ))).^0.25;
    etaco  = (eta(icz(1:end-1),icx(1:end-1)).*eta(icz(2:end),icx(1:end-1)) ...
           .* eta(icz(1:end-1),icx(2:end  )).*eta(icz(2:end),icx(2:end  ))).^0.25;
    tupdco = (tupd(icz(1:end-1),icx(1:end-1)).*tupd(icz(2:end),icx(1:end-1)) ...
           .* tupd(icz(1:end-1),icx(2:end  )).*tupd(icz(2:end),icx(2:end  ))).^0.25;
elseif meansw == 1 % Arithmetic
    etavpco  = (etavp(icz(1:end-1),icx(1:end-1))+etavp(icz(2:end),icx(1:end-1)) ...
             +  etavp(icz(1:end-1),icx(2:end  ))+etavp(icz(2:end),icx(2:end  ))).*0.25;
    etaco  = (eta(icz(1:end-1),icx(1:end-1))+eta(icz(2:end),icx(1:end-1)) ...
           +  eta(icz(1:end-1),icx(2:end  ))+eta(icz(2:end),icx(2:end  ))).*0.25;
    tupdco = (tupd(icz(1:end-1),icx(1:end-1))+tupd(icz(2:end),icx(1:end-1)) ...
           +  tupd(icz(1:end-1),icx(2:end  ))+tupd(icz(2:end),icx(2:end  ))).*0.25;
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
deltac = sqrt(KD.*min(1e26,zetavp));

% heat dissipation (entropy production) rate
if Nz==1 && Nx==1
    diss = 0.*T;  % no dissipation in 0-D mode (no diffusion, no shear deformation, no segregation)
else
    [grdTx ,grdTz ] = gradient(T(icz,icx),h);
    exz_ce = (exz(1:end-1,1:end-1)+exz(2:end,1:end-1)+exz(1:end-1,2:end)+exz(2:end,2:end))./4;
    diss = min(1e-3,diss ...
         + kT./T.*(grdTz (2:end-1,2:end-1).^2 + grdTx (2:end-1,2:end-1).^2) ...
         + eta.*exx.^2 + eta.*ezz.^2 + 2.*eta.*exz_ce.^2 ...
         + (mu+1e-9)./Ks .* ((wm(1:end-1,2:end-1)+wm(2:end,2:end-1))./2).^2 ...
         + (mu+1e-9)./Ks .* ((um(2:end-1,1:end-1)+um(2:end-1,2:end))./2).^2)/2;
end

% update time step
dtk = (h/2)^2/max([kT(:)./rho(:)./cP(:);kd(:)]);                           % diffusive time step size
dta =  h/2   /max(abs([Um(:);Wm(:);Ux(:);Wx(:)]));                         % advective time step size
dt  = min([1.5*dto,min(CFL*[dtk,dta]),dtmax]);     % time step size

UDtime = UDtime + toc;
