%*** update phase equilibrium
eqtime = tic;

var.c      = reshape(c,Nx*Nz,cal.ncmp);   % component fractions [wt]
var.T      = reshape(T,Nx*Nz,1)-273.15;   % temperature [C]
var.P      = reshape(Pt,Nx*Nz,1)/1e9;     % pressure [GPa]
var.m      = reshape(mq,Nx*Nz,1);         % melt fraction [wt]
var.H2O    = var.c(:,end);                % water concentration [wt]
var.X      = reshape(cm_oxd_all,Nz*Nx,9); % melt oxide fractions [wt %]
cal.H2Osat = fluidsat(var);               % water saturation [wt]

[var,cal]  = leappart(var,cal,'E');

Tsol   = reshape(cal.Tsol,Nz,Nx);
Tliq   = reshape(cal.Tliq,Nz,Nx);
H2Osat = reshape(cal.H2Osat,Nz,Nx);

mq = reshape(var.m.*(var.m>eps^0.5),Nz,Nx);
xq = reshape(var.x.*(var.x>eps^0.5),Nz,Nx);
mq = mq./(mq+xq);
xq = xq./(mq+xq);

cxq = reshape(var.cx,Nz,Nx,cal.ncmp); cxq = cxq./sum(cxq,3);
cmq = reshape(var.cm,Nz,Nx,cal.ncmp); cmq = cmq./sum(cmq,3);

Kx  = reshape(cal.Kx,Nz,Nx,cal.ncmp);

Ktrc  = zeros(Nz,Nx,cal.ntrc);
for i = 1:cal.ntrc  
    % update bulk partitioning coefficients
    for j=1:cal.nmem; Ktrc(:,:,i) = Ktrc(:,:,i) + cal.Ktrc_mem(i,j) .* cx_mem(:,:,j)./100; end
end  
% update trace element phase compositions
trcmq = trc./(m + x.*Ktrc);
trcxq = trc./(m./Ktrc + x);

% disequilibrium - equilibrium transition at low melt fraction
f_diseq = 1 - exp(-mu./mutra) - exp(-chi./mutra);

% equilibrium component transfer rates
Gmc_eq = (cmq.*mq.*rho-Cm)/(tau_phs+3*dt);
Gxc_eq = -Gmc_eq;

Gmtrc_eq = (trcmq.*mq.*rho-TRCm)/(tau_phs+3*dt);
Gxtrc_eq = -Gmtrc_eq;

% net phase mass transfer rates
Gm = sum(Gmc_eq,3);
Gx = -Gm;

% disequlibrium apporach
% 1. fractional melting/crystallisation compositions
cxqf = cm.*Kx; cxqf = cxqf./sum(cxqf,3);
cmqf = cx./Kx; cmqf = cmqf./sum(cmqf,3);

trcxqf = trcm.*Ktrc;
trcmqf = trcx./Ktrc;

% 2. phase change component transfer (sum over comp to Gm, Gx, sum over phs to 0)
Gmc_phs = max(0,Gm).*cmqf + min(0,Gm).*cxqf;
Gxc_phs = -Gmc_phs;

Gmtrc_phs = max(0,Gm).*trcmqf + min(0,Gm).*trcxqf;
Gxtrc_phs = -Gmtrc_phs;

% 3. chemical exchange component transfer (sum over comp to 0, sum over phs to 0)
Gmc_cex = Gmc_eq - Gmc_phs;
Gxc_cex = Gxc_eq - Gxc_phs;

Gmtrc_cex = Gmtrc_eq - Gmtrc_phs;
Gxtrc_cex = Gxtrc_eq - Gxtrc_phs;

% 4. combined disequilibrium component tranfer rates
Gmc_diseq = (Gmc_phs + Gmc_cex);
Gxc_diseq = (Gxc_phs + Gxc_cex);

Gmtrc_diseq = (Gmtrc_phs + Gmtrc_cex);
Gxtrc_diseq = (Gxtrc_phs + Gxtrc_cex);

% 5. Adding equlibrum and disequilibrium components weighted by smooth
% transition functon
Gmc = (1-f_diseq).*Gmc_eq + f_diseq .* Gmc_diseq;
Gxc = -Gmc;

Gmtrc = (1-f_diseq).*Gmtrc_eq + f_diseq .* Gmtrc_diseq;
Gxtrc = -Gmtrc;


% extract, extrude, and intrude melt

findmoho;

% Gem  = min(0,mthr-m).*RHO/(tau_e+3*dt).*(ZZ<repmat(LAB_depth+2e3,Nz,1));%(etax>1e20);
Gem  = min(0,mthr-m).*RHO/(tau_e+3*dt).*(etax>1e20);
Gemc = cm  .*Gem;
Gemt = trcm.*Gem;
Gems = sm  .*Gem;

findmelt;

%Extrusion

extr_shape = (1-path_ratio).*exp(-abs(ZZ - h/2) / bnd_w);
path_shape = path_ratio.*(ZZ<max(h,melt_depth));

Gex = erupt_ratio * (extr_shape+path_shape).*(-sum(Gem,1))./sum((extr_shape+path_shape),1);
% for i=1:2; Gex = Gex + diff(Gex(:,icx),2,2)./4; end

Gexc = erupt_ratio * (extr_shape+path_shape).*(-sum(Gemc,1))./sum((extr_shape+path_shape),1);
% for i=1:2; Gexc = Gexc + diff(Gexc(:,icx,:),2,2)./4; end

Gext = erupt_ratio * (extr_shape+path_shape).*(-sum(Gemt,1))./sum((extr_shape+path_shape),1);
% for i=1:2; Gext = Gext + diff(Gext(:,icx,:),2,2)./4; end

sx0   = cPx.*log((T0+273.15)./Tref);
Gems0 = sx0 .* Gem;
Gexs  = erupt_ratio * extr_shape.*(-sum(Gems0,1))./sum(extr_shape,1) ...
      + erupt_ratio * path_shape.*(-sum(Gems ,1))./sum(path_shape,1);
% for i=1:2; Gexs = Gexs + diff(Gexs(:,icx,:),2,2)./4; end

% Intrusion

intr_shape = (1-path_ratio).*exp(-abs(ZZ - min(MOHO_depth,melt_depth-h)) / bnd_w) + path_ratio.*(ZZ>min(MOHO_depth,melt_depth-h) & ZZ<melt_depth);
Gin = (1-erupt_ratio) * intr_shape.*(-sum(Gem,1))./sum(intr_shape,1);
% for i=1:2; Gin = Gin + diff(Gin(:,icx),2,2)./4; end
Gex = Gex + Gin.*(T-273.15<=Tsol);
Gin =       Gin.*(T-273.15> Tsol);

Ginc = (1-erupt_ratio) * intr_shape.*(-sum(Gemc,1))./sum(intr_shape,1);
% for i=1:2; Ginc = Ginc + diff(Ginc(:,icx,:),2,2)./4; end
Gexc = Gexc + Ginc.*(T-273.15<=Tsol);
Ginc =        Ginc.*(T-273.15> Tsol);

Gint = (1-erupt_ratio) * intr_shape.*(-sum(Gemt,1))./sum(intr_shape,1);
% for i=1:2; Gint = Gint + diff(Gint(:,icx,:),2,2)./4; end
Gext = Gext + Gint.*(T-273.15<=Tsol);
Gint =        Gint.*(T-273.15> Tsol);

Gins = (1-erupt_ratio) * intr_shape.*(-sum(Gems,1))./sum(intr_shape,1);
% for i=1:2; Gins = Gins + diff(Gins(:,icx,:),2,2)./4; end
Gexs = Gexs + Gins.*(T-273.15<=Tsol);
Gins =        Gins.*(T-273.15> Tsol);

eqtime = toc(eqtime);
EQtime = EQtime + eqtime;