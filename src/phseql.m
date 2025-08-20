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

% phase mass transfer rates 
Gm  = (mq-m).*RHO/(tau_r+3*dt);
Gx  = (xq-x).*RHO/(tau_r+3*dt); 

Gmc = (cmq.*mq-cm.*m).*RHO/(tau_r+3*dt);
Gxc = (cxq.*xq-cx.*x).*RHO/(tau_r+3*dt);

% extract and erupt melt

findmoho

 Gem = min(0,mthr-m).*RHO/(tau_e+3*dt);

% Eruption

 Gex_erupt = topshape.*(-sum(Gem,1))./sum(topshape,1);
 for i=1:4; Gex_erupt = Gex_erupt + diff(Gex_erupt(:,icx),2,2)./8; end

 Gemc = cm.*Gem;
 Gexc_erupt = topshape.*(-sum(Gemc,1))./sum(topshape,1);
 for i=1:4; Gexc_erupt = Gexc_erupt + diff(Gexc_erupt(:,icx,:),2,2)./8; end

 Gemt = trcm.*Gem;
 Gext_erupt = topshape.*(-sum(Gemt,1))./sum(topshape,1);
 for i=1:4; Gext_erupt = Gext_erupt + diff(Gext_erupt(:,icx,:),2,2)./8; end

% Emplacement 

 mohoshape = exp( -abs(ZZ - moho_depth) / bnd_w );
 Gex_intr = mohoshape.*(-sum(Gem,1))./sum(mohoshape,1);
 for i=1:4; Gex_intr = Gex_intr + diff(Gex_intr(:,icx),2,2)./8; end

 Gemc = cm.*Gem;
 Gexc_intr = mohoshape.*(-sum(Gemc,1))./sum(mohoshape,1);
 for i=1:4; Gexc_intr = Gexc_intr + diff(Gexc_intr(:,icx,:),2,2)./8; end

 Gemt = trcm.*Gem;
 Gext_intr = mohoshape.*(-sum(Gemt,1))./sum(mohoshape,1);
 for i=1:4; Gext_intr = Gext_intr + diff(Gext_intr(:,icx,:),2,2)./8; end


%Combining the Eruption and Emplacement sections based on the ratio

Gex  = erupt_ratio * Gex_erupt  + (1-erupt_ratio) * Gex_intr;
Gexc = erupt_ratio * Gexc_erupt + (1-erupt_ratio) * Gexc_intr;
Gext = erupt_ratio * Gext_erupt + (1-erupt_ratio) * Gext_intr;

cxq = reshape(var.cx,Nz,Nx,cal.ncmp);
cmq = reshape(var.cm,Nz,Nx,cal.ncmp);

eqtime = toc(eqtime);
EQtime = EQtime + eqtime;

