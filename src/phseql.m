%*** update phase equilibrium
eqtime = tic;

Ktrc  = zeros(Nz,Nx,cal.ntrc);
for i = 1:cal.ntrc  
    % update bulk partitioning coefficients
    for j=1:cal.nmem; Ktrc(:,:,i) = Ktrc(:,:,i) + cal.Ktrc_mem(i,j) .* majx_mem(:,:,j)./100; end
end
cal.Ktrc = reshape(Ktrc,Nz*Nx,cal.ntrc);

var.c      = reshape(c,Nx*Nz,cal.ncmp+cal.ntrc);   % component fractions [wt]
var.T      = reshape(T,Nx*Nz,1)-273.15;   % temperature [C]
var.P      = reshape(Pt,Nx*Nz,1)/1e9;     % pressure [GPa]
var.m      = reshape(mq,Nx*Nz,1);         % melt fraction [wt]
var.X      = reshape(majm_oxd_all,Nz*Nx,9); % melt oxide fractions [wt %]
cal.H2Osat = fluidsat(var);               % water saturation [wt]

[var,cal]  = leappart(var,cal,'E');

Tsol   = reshape(cal.Tsol,Nz,Nx);
Tliq   = reshape(cal.Tliq,Nz,Nx);
H2Osat = reshape(cal.H2Osat,Nz,Nx);

mq = reshape(var.m.*(var.m>eps^0.5),Nz,Nx);
xq = reshape(var.x.*(var.x>eps^0.5),Nz,Nx);
mq = mq./(mq+xq);
xq = xq./(mq+xq);
x  = xq;  m = mq;

cxq = reshape(var.cx,Nz,Nx,cal.ncmp+cal.ntrc); cxq = cxq./sum(cxq,3);
cmq = reshape(var.cm,Nz,Nx,cal.ncmp+cal.ntrc); cmq = cmq./sum(cmq,3);
majmq = cmq(:,:,1:cal.ncmp);  majmq = majmq./sum(majmq,3);
majxq = cxq(:,:,1:cal.ncmp);  majxq = majxq./sum(majxq,3);
trcmq = cmq(:,:,cal.ncmp+1:end)/trcscale;
trcxq = cxq(:,:,cal.ncmp+1:end)/trcscale;

Kx = reshape(cal.Kx,Nz,Nx,cal.ncmp+cal.ntrc);

% disequilibrium - equilibrium transition at low melt fraction
f_diseq = 1 - exp(-mu./mutra) - exp(-chi./mutra); %(1+erf((mu-mutra)./mutra*4))/2;%

% equilibrium component transfer rates
Gmc_eq  = (cmq.*mq.*RHO-Cm)/(tau_phs+3*dt);

% net phase mass transfer rates
Gm = sum(Gmc_eq,3);
Gx = -Gm;

% disequlibrium apporach
% 1. fractional melting/crystallisation compositions
cxqf = cm.*Kx; cxqf = cxqf./sum(cxqf,3);
cmqf = cx./Kx; cmqf = cmqf./sum(cmqf,3);

% 2. phase change component transfer (sum over comp to Gm, Gx, sum over phs to 0)
Gmc_phs = max(0,Gm).*cmqf   + min(0,Gm).*cxqf  ;

% 3. chemical exchange component transfer (sum over comp to 0, sum over phs to 0)
Gmc_cex   = Gmc_eq   - Gmc_phs  ;

% 4. combined disequilibrium component tranfer rates
Gmc_diseq = (Gmc_phs   + Gmc_cex);

% 5. Adding equlibrum and disequilibrium components weighted by smooth
% transition functon
Gmc = (1-f_diseq).*Gmc_eq + f_diseq .* Gmc_diseq;
Gxc = -Gmc;

% extract, extrude, and intrude melt

findmoho;

% Gem  = min(0,mthr-m).*RHO/(tau_e+3*dt).*(ZZ<repmat(LAB_depth+2e3,Nz,1));%(etax>1e20);
Gem  = min(0,mthr    .*RHO-M )/(tau_e+3*dt).*(etax>1e20);
Gemc = min(0,mthr.*cm.*RHO-Cm)/(tau_e+3*dt).*(etax>1e20);
Gems = sm.*Gem;

findmelt;

%Extrusion

extr_shape = (1-path_ratio).*exp(-abs(ZZ - h/2) / bnd_w);
path_shape = path_ratio.*(ZZ<max(h,melt_depth));

Gex  = erupt_ratio * (extr_shape+path_shape).*(-sum(Gem ,1))./sum((extr_shape+path_shape),1);

Gexc = erupt_ratio * (extr_shape+path_shape).*(-sum(Gemc,1))./sum((extr_shape+path_shape),1);

sx0   = cPx.*log((T0+273.15)./Tref);
Gems0 = sx0 .* Gem;
Gexs  = erupt_ratio * extr_shape.*(-sum(Gems0,1))./sum(extr_shape,1) ...
      + erupt_ratio * path_shape.*(-sum(Gems ,1))./sum(path_shape,1);

% Intrusion

intr_shape = (1-path_ratio).*exp(-abs(ZZ - min(MOHO_depth,melt_depth-h)) / bnd_w) + path_ratio.*(ZZ>min(MOHO_depth,melt_depth-h) & ZZ<melt_depth);
Gin = (1-erupt_ratio) * intr_shape.*(-sum(Gem,1))./sum(intr_shape,1);
Gex = Gex + Gin.*(T-273.15<=Tsol);
Gin =       Gin.*(T-273.15> Tsol);

Ginc = (1-erupt_ratio) * intr_shape.*(-sum(Gemc,1))./sum(intr_shape,1);
Gexc = Gexc + Ginc.*(T-273.15<=Tsol);
Ginc =        Ginc.*(T-273.15> Tsol);

Gins = (1-erupt_ratio) * intr_shape.*(-sum(Gems,1))./sum(intr_shape,1);
Gexs = Gexs + Gins.*(T-273.15<=Tsol);
Gins =        Gins.*(T-273.15> Tsol);

eqtime = toc(eqtime);
EQtime = EQtime + eqtime;