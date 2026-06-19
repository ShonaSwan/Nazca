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

cxq = reshape(var.cx,Nz,Nx,cal.ncmp);
cmq = reshape(var.cm,Nz,Nx,cal.ncmp);

Ktrc    = zeros(Nz,Nx,cal.ntrc);
trcmq   = zeros(Nz,Nx,cal.ntrc);
trcxq   = zeros(Nz,Nx,cal.ntrc);

for i = 1:cal.ntrc  
    % update bulk partitioning coefficients
    for j=1:cal.nmem; Ktrc(:,:,i) = Ktrc(:,:,i) + cal.Ktrc_mem(i,j) .* cx_mem(:,:,j)./100; end

    % update trace element phase compositions
    trcmq(:,:,i) = trc(:,:,i)./(m + x.*Ktrc(:,:,i));
    trcxq(:,:,i) = trc(:,:,i)./(m./Ktrc(:,:,i) + x);
end  

% smooth weighting for small melt fraction (Low melt → force equilibrium)
melt_weight = 0.5 * (1 + tanh((m - 0.005)/0.003)); 

% phase mass transfer rates
Gm  = (mq-m).*RHO/(tau_phs+3*dt);
Gx  = (xq-x).*RHO/(tau_phs+3*dt);

% equilibrium component transfer rates
Gmcq = (cmq.*mq-cm.*m).*RHO/(tau_phs+3*dt);
Gxcq = (cxq.*xq-cx.*x).*RHO/(tau_phs+3*dt);

% disequilibrium component transfer rates

% phase change component transfer (sum over comp to Gm, Gx, sum over phs to 0)
Gmc = max(0,Gm).*cmq - max(0,Gx).*cxq;
Gxc = max(0,Gx).*cxq - max(0,Gm).*cmq;

Gmcex = (cmq - cm) .* M / tau_cex;
Gxcex = (cxq - cx) .* X / tau_cex;

% chemical exchange component transfer (sum over comp to 0, sum over phs to 0)
Gmtrc = max(0,Gm).*trcmq - max(0,Gx).*trcxq;
Gxtrc = max(0,Gx).*trcxq - max(0,Gm).*trcmq;

Gmtrcex = (trcmq - trcm) .* M / tau_cex;
Gxtrcex = (trcxq - trcx) .* X / tau_cex;


% extract, extrude, and intrude melt

findmoho;

Gem  = min(0,mthr-m).*RHO/(tau_e+3*dt).*(ZZ<repmat(LAB_depth+2e3,Nz,1));%(etax>1e20);
Gemc = cm  .*Gem;
Gemt = trcm.*Gem;
Gems = sm  .*Gem;

findmelt;

%Extrusion

extr_shape = (1-path_ratio).*exp(-abs(ZZ - h/2) / bnd_w);
path_shape = path_ratio.*(ZZ<max(h,melt_depth));
Gex = erupt_ratio * (extr_shape+path_shape).*(-sum(Gem,1))./sum((extr_shape+path_shape),1);
for i=1:2; Gex = Gex + diff(Gex(:,icx),2,2)./4; end

Gexc = erupt_ratio * (extr_shape+path_shape).*(-sum(Gemc,1))./sum((extr_shape+path_shape),1);
for i=1:2; Gexc = Gexc + diff(Gexc(:,icx,:),2,2)./4; end

Gext = erupt_ratio * (extr_shape+path_shape).*(-sum(Gemt,1))./sum((extr_shape+path_shape),1);
for i=1:2; Gext = Gext + diff(Gext(:,icx,:),2,2)./4; end

sx0   = cPx.*log((T0+273.15)./Tref);
Gems0 = sx0 .* Gem;
Gexs  = erupt_ratio * extr_shape.*(-sum(Gems0,1))./sum(extr_shape,1) ...
      + erupt_ratio * path_shape.*(-sum(Gems ,1))./sum(path_shape,1);
for i=1:2; Gexs = Gexs + diff(Gexs(:,icx,:),2,2)./4; end

% Intrusion

intr_shape = (1-path_ratio).*exp(-abs(ZZ - min(MOHO_depth,melt_depth-h)) / bnd_w) + path_ratio.*(ZZ>min(MOHO_depth,melt_depth-h) & ZZ<melt_depth);
Gin = (1-erupt_ratio) * intr_shape.*(-sum(Gem,1))./sum(intr_shape,1);
for i=1:2; Gin = Gin + diff(Gin(:,icx),2,2)./4; end
Gex = Gex + Gin.*(T-273.15<=Tsol);
Gin =       Gin.*(T-273.15> Tsol);

Ginc = (1-erupt_ratio) * intr_shape.*(-sum(Gemc,1))./sum(intr_shape,1);
for i=1:2; Ginc = Ginc + diff(Ginc(:,icx,:),2,2)./4; end
Gexc = Gexc + Ginc.*(T-273.15<=Tsol);
Ginc =        Ginc.*(T-273.15> Tsol);

Gint = (1-erupt_ratio) * intr_shape.*(-sum(Gemt,1))./sum(intr_shape,1);
for i=1:2; Gint = Gint + diff(Gint(:,icx,:),2,2)./4; end
Gext = Gext + Gint.*(T-273.15<=Tsol);
Gint =        Gint.*(T-273.15> Tsol);

Gins = (1-erupt_ratio) * intr_shape.*(-sum(Gems,1))./sum(intr_shape,1);
for i=1:2; Gins = Gins + diff(Gins(:,icx,:),2,2)./4; end
Gexs = Gexs + Gins.*(T-273.15<=Tsol);
Gins =        Gins.*(T-273.15> Tsol);

eqtime = toc(eqtime);
EQtime = EQtime + eqtime;