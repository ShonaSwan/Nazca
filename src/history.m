% record run history

dsumMdtoo = dsumMdto; dsumMdto = dsumMdt;
dsumSdtoo = dsumSdto; dsumSdto = dsumSdt;
dsumCdtoo = dsumCdto; dsumCdto = dsumCdt;

stp = max(1,step);

% record model time
hist.time(stp) = time;

% record total mass, heat, component mass in model (assume hy = 1, unit length in third dimension)
hist.sumM(stp  ) = sum(rho(:)*h*h*1);  % [kg]
hist.sumS(stp  ) = sum(  S(:)*h*h*1);  % [J]
hist.sumC(stp,:) = squeeze(sum(sum(C*h*h*1,1),2))+eps; % [kg]

% record expected rates of change by volume change and imposed boundaries layers
dsumMdt = sum(X(1,:).*Wx(1,2:end-1)*h*1) - sum(X(end,:).*Wx(end,2:end-1)*h*1) ...
        + sum(X(:,1).*Ux(2:end-1,1)*h*1) - sum(X(:,end).*Ux(2:end-1,end)*h*1) ...
        + sum(M(1,:).*Wm(1,2:end-1)*h*1) - sum(M(end,:).*Wm(end,2:end-1)*h*1) ...
        + sum(M(:,1).*Um(2:end-1,1)*h*1) - sum(M(:,end).*Um(2:end-1,end)*h*1);  % [kg/s]
dsumSdt = sum(sum(bnd_S*h*h*1)) + sum(sum(diss_h*h*h*1)) ...
        + sum(X(1,:).*sx(1,:).*Wx(1,2:end-1)*h*1) - sum(X(end,:).*sx(end,:).*Wx(end,2:end-1)*h*1) ...
        + sum(X(:,1).*sx(:,1).*Ux(2:end-1,1)*h*1) - sum(X(:,end).*sx(:,end).*Ux(2:end-1,end)*h*1) ...
        + sum(M(1,:).*sm(1,:).*Wm(1,2:end-1)*h*1) - sum(M(end,:).*sm(end,:).*Wm(end,2:end-1)*h*1) ...
        + sum(M(:,1).*sm(:,1).*Um(2:end-1,1)*h*1) - sum(M(:,end).*sm(:,end).*Um(2:end-1,end)*h*1);  % [J /s]
dsumCdt = squeeze(sum(sum(bnd_C*h*h*1,1),2) ...
        + sum(X(1,:).*cx(1,:,:).*Wx(1,2:end-1)*h*1,2) - sum(X(end,:).*cx(end,:,:).*Wx(end,2:end-1)*h*1,2) ...
        + sum(X(:,1).*cx(:,1,:).*Ux(2:end-1,1)*h*1,1) - sum(X(:,end).*cx(:,end,:).*Ux(2:end-1,end)*h*1,1) ...
        + sum(M(1,:).*cm(1,:,:).*Wm(1,2:end-1)*h*1,2) - sum(M(end,:).*cm(end,:,:).*Wm(end,2:end-1)*h*1,2) ...
        + sum(M(:,1).*cm(:,1,:).*Um(2:end-1,1)*h*1,1) - sum(M(:,end).*cm(:,end,:).*Um(2:end-1,end)*h*1,1)).';  % [kg/s]

if step>=2; hist.DM(stp  ) = (a2*hist.DM(max(1,stp-1)  ) + a3*hist.DM(max(1,stp-2)  ) + (b1*dsumMdt + b2*dsumMdto + b3*dsumMdtoo)*dt)/a1; else; hist.DM(stp  ) = 0; end  % [kg]
if step>=2; hist.DS(stp  ) = (a2*hist.DS(max(1,stp-1)  ) + a3*hist.DS(max(1,stp-2)  ) + (b1*dsumSdt + b2*dsumSdto + b3*dsumSdtoo)*dt)/a1; else; hist.DS(stp  ) = 0; end  % [kg]
if step>=2; hist.DC(stp,:) = (a2*hist.DC(max(1,stp-1),:) + a3*hist.DC(max(1,stp-2),:) + (b1*dsumCdt + b2*dsumCdto + b3*dsumCdtoo)*dt)/a1; else; hist.DC(stp,:) = zeros(1,cal.ncmp); end  % [kg]

% if stp>=2; hist.DM(stp  ) = hist.DM(max(1,stp-1)  ) + dsumMdt*dt; else; hist.DM(stp  ) = 0; end  % [kg]
% if stp>=2; hist.DS(stp  ) = hist.DS(max(1,stp-1)  ) + dsumSdt*dt; else; hist.DS(stp  ) = 0; end  % [kg]
% if stp>=2; hist.DC(stp,:) = hist.DC(max(1,stp-1),:) + dsumCdt*dt; else; hist.DC(stp,:) = zeros(1,cal.ncmp); end  % [kg]

% record conservation error of mass M, heat S, components C
hist.EM(stp  ) = (hist.sumM(stp  ) - hist.DM(stp  ))./hist.sumM(1  ) - 1;  % [kg/kg]
hist.ES(stp  ) = (hist.sumS(stp  ) - hist.DS(stp  ))./hist.sumS(1  ) - 1;  % [JK/JK]
hist.EC(stp,:) = (hist.sumC(stp,:) - hist.DC(stp,:))./hist.sumC(1,:) - 1;  % [kg/kg]

% if step==2; hist.sumM(1) = hist.sumM(2); end

% record variable and coefficient diagnostics
hist.W(stp,1) = min(min(-W(:,2:end-1)));
hist.W(stp,2) = mean(mean(abs(W(:,2:end-1))));
hist.W(stp,3) = max(max(-W(:,2:end-1)));

hist.U(stp,1) = min(min(U(2:end-1,:)));
hist.U(stp,2) = mean(mean(abs(U(2:end-1,:))));
hist.U(stp,3) = max(max(U(2:end-1,:)));

hist.Pf(stp,1) = min(min(Pf(2:end-1,2:end-1)));
hist.Pf(stp,2) = mean(mean(abs(Pf(2:end-1,2:end-1))));
hist.Pf(stp,3) = max(max(Pf(2:end-1,2:end-1)));

hist.Pchmb(stp,1) = Pchmb;

hist.x(stp,1) = min(min(x));
hist.x(stp,2) = mean(mean(x));
hist.x(stp,3) = max(max(x));

hist.m(stp,1) = min(min(m));
hist.m(stp,2) = mean(mean(m));
hist.m(stp,3) = max(max(m));

hist.chi(stp,1) = min(min(chi));
hist.chi(stp,2) = mean(mean(chi));
hist.chi(stp,3) = max(max(chi));

hist.mu(stp,1) = min(min(mu));
hist.mu(stp,2) = mean(mean(mu));
hist.mu(stp,3) = max(max(mu));

hist.T(stp,1) = min(min(T));
hist.T(stp,2) = mean(mean(T));
hist.T(stp,3) = max(max(T));

hist.Tsol(stp,1) = min(min(Tsol));
hist.Tsol(stp,2) = mean(mean(Tsol));
hist.Tsol(stp,3) = max(max(Tsol));

hist.Tliq(stp,1) = min(min(Tliq));
hist.Tliq(stp,2) = mean(mean(Tliq));
hist.Tliq(stp,3) = max(max(Tliq));

for i=1:cal.ncmp
    hist.c(stp,1,i) = min(min(c(:,:,i)));
    hist.c(stp,2,i) = mean(mean(c(:,:,i)));
    hist.c(stp,3,i) = max(max(c(:,:,i)));
end
for i=1:cal.noxd
    hist.c_oxd(stp,1,i) = min(min(c_oxd(:,:,i)));
    hist.c_oxd(stp,2,i) = mean(mean(c_oxd(:,:,i)));
    hist.c_oxd(stp,3,i) = max(max(c_oxd(:,:,i)));
end

for i=1:cal.ncmp
    hist.cx(stp,1,i) = min(min(cx(:,:,i)));
    hist.cx(stp,2,i) = sum(sum(cx(:,:,i).*x.*rho))./sum(sum(x.*rho));
    hist.cx(stp,3,i) = max(max(cx(:,:,i)));
end
for i=1:cal.noxd
    hist.cx_oxd(stp,1,i) = min(min(cx_oxd(:,:,i)));
    hist.cx_oxd(stp,2,i) = sum(sum(cx_oxd(:,:,i).*x.*rho))./sum(sum(x.*rho));
    hist.cx_oxd(stp,3,i) = max(max(cx_oxd(:,:,i)));
end
for i=1:cal.nmsy
    hist.cx_msy(stp,1,i) = min(min(cx_msy(:,:,i)));
    hist.cx_msy(stp,2,i) = sum(sum(cx_msy(:,:,i).*x.*rho))./sum(sum(x.*rho));
    hist.cx_msy(stp,3,i) = max(max(cx_msy(:,:,i)));
end

hist.rhox(stp,1) = min(min(rhox));
hist.rhox(stp,2) = sum(sum(rhox.*x.*rho))./sum(sum(x.*rho));
hist.rhox(stp,3) = max(max(rhox));

for i=1:cal.ncmp
    hist.cm(stp,1,i) = min(min(cm(:,:,i)));
    hist.cm(stp,2,i) = sum(sum(cm(:,:,i).*m.*rho))./sum(sum(m.*rho));
    hist.cm(stp,3,i) = max(max(cm(:,:,i)));
end
for i=1:cal.noxd
    hist.cm_oxd(stp,1,i) = min(min(cm_oxd(:,:,i)));
    hist.cm_oxd(stp,2,i) = sum(sum(cm_oxd(:,:,i).*m.*rho))./sum(sum(m.*rho));
    hist.cm_oxd(stp,3,i) = max(max(cm_oxd(:,:,i)));
end

hist.rhom(stp,1) = min(min(rhom));
hist.rhom(stp,2) = sum(sum(rhom.*m))./sum(sum(m));
hist.rhom(stp,3) = max(max(rhom));

hist.etam(stp,1) = min(min(etam));
hist.etam(stp,2) = sum(sum(etam.*m))./sum(sum(m));
hist.etam(stp,3) = max(max(etam));

hist.Gm(stp,1) = min(min(Gm));
hist.Gm(stp,2) = mean(mean(Gm));
hist.Gm(stp,3) = max(max(Gm));

hist.Gx(stp,1) = min(min(Gx));
hist.Gx(stp,2) = mean(mean(Gx));
hist.Gx(stp,3) = max(max(Gx));

hist.dV(stp,1) = min(min(Div_V));
hist.dV(stp,2) = mean(mean(Div_V));
hist.dV(stp,3) = max(max(Div_V));

hist.rho(stp,1) = min(min(rho));
hist.rho(stp,2) = mean(mean(rho));
hist.rho(stp,3) = max(max(rho));

hist.eta(stp,1) = min(min(eta));
hist.eta(stp,2) = geomean(geomean(eta));
hist.eta(stp,3) = max(max(eta));

hist.wx(stp,1) = min(min(-(chi([1,1:end],:)+chi([1:end,end],:))/2.*wx(:,2:end-1)));
hist.wx(stp,2) = mean(mean(abs((chi([1,1:end],:)+chi([1:end,end],:))/2.*wx(:,2:end-1))));
hist.wx(stp,3) = max(max(-(chi([1,1:end],:)+chi([1:end,end],:))/2.*wx(:,2:end-1)));

hist.wm(stp,1) = min(min(-(mu([1,1:end],:)+mu([1:end,end],:))/2.*wm(:,2:end-1)));
hist.wm(stp,2) = mean(mean(abs((mu([1,1:end],:)+mu([1:end,end],:))/2.*wm(:,2:end-1))));
hist.wm(stp,3) = max(max(-(mu([1,1:end],:)+mu([1:end,end],:))/2.*wm(:,2:end-1)));

hist.dH(stp,1) = -sum(bnd_S(1:round(Nz/2)  ,:).*T(1:round(Nz/2)  ,:).*h^2,'all')/L;
hist.dH(stp,2) = -sum(bnd_S(round(Nz/2):end,:).*T(round(Nz/2):end,:).*h^2,'all')/L;
hist.dH(stp,3) = -sum(bnd_S.*T.*h^2,'all')/L;