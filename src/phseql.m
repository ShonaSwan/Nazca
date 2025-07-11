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
Gm  = (mq-m).*RHO/max(tau_r,10*dt);
Gx  = (xq-x).*RHO/max(tau_r,10*dt); 

Gmc = (cmq.*mq-cm.*m).*RHO/max(tau_r,10*dt);
Gxc = (cxq.*xq-cx.*x).*RHO/max(tau_r,10*dt);

% extract and erupt melt
Gem = min(0,mthr-m).*RHO/max(tau_e,10*dt);
Gex = topshape.*(-sum(Gem,1))./sum(topshape,1);
for i=1:4; Gex = Gex + diffus(Gex,1/8*ones(size(rp)),1,[1,2],BCD); end

Gemc = cm.*Gem;
Gexc = topshape.*(-sum(Gemc,1))./sum(topshape,1);
for i=1:5; Gexc = Gexc + diffus(Gexc,1/8*ones(size(rp)),1,[1,2],BCD); end

Gemt = trcm.*Gem;
Gext = topshape.*(-sum(Gemt,1))./sum(topshape,1);
for i=1:4; Gext = Gext + diffus(Gext,1/8*ones(size(rp)),1,[1,2],BCD); end

cxq = reshape(var.cx,Nz,Nx,cal.ncmp);
cmq = reshape(var.cm,Nz,Nx,cal.ncmp);

eqtime = toc(eqtime);
EQtime = EQtime + eqtime;