function [datafit,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmpfit,MLT_cmpfit,SYS_cmpfit] = ModelFitNP(model,T,P,SYS,PHS,cal,reg)

% get number of data points
np = length(T);

% get model parameters from model vector
cal.T0    = model(               (1:cal.ncmp-1)).';
cal.A     = model(1*(cal.ncmp-1)+(1:cal.ncmp-1)).';
cal.B     = model(2*(cal.ncmp-1)+(1:cal.ncmp-1)).';
cal.r     = model(3*(cal.ncmp-1)+(1:cal.ncmp-1)).';
cal.dTH2O = 1400 * 1200./cal.T0; %model(4*(cal.ncmp-1)+(1:cal.ncmp-1)).';
cmp_mem   = reshape(model(5*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem)),cal.ncmp,cal.nmem);
cmp_oxd   = cmp_mem*cal.mem_oxd./100;

% find system component composition by least-squares
SYS_cmpfit = lsqregcmp(cmp_oxd(1:end-1,1:end-1).',SYS(:,1:end-1),reg);
SYS_cmpfit = [SYS_cmpfit,SYS(:,end)/100];
SYS_cmpfit = SYS_cmpfit./sum(SYS_cmpfit,2);

% get fitted phase oxide compositions
SYS_oxdfit = SYS_cmpfit*cmp_oxd;

% get local phase equilibrium
var.m = ones(np,1)/2; var.x = var.m; var.f = 0*var.m;
var.c      = SYS_cmpfit;        % component fractions [wt]
% var.X      = var.c*cmp_oxd; % melt oxide fractions [wt %]
var.T      = T;             % temperature [C]
var.P      = P/10;         % pressure [GPa]
var.H2O    = SYS_cmpfit(:,end); % water concentration [wt]
cal.H2Osat = var.H2O./(PHS(:,1)/100);%fluidsat(var);
[var,cal]  = leappart(var,cal,'E');

% get fitted phase oxide compositions
MLT_cmpfit = var.cm;
SOL_cmpfit = var.cx;
MLT_oxdfit = var.cm*cmp_oxd;
SOL_oxdfit = var.cx*cmp_oxd;
SOL_memfit = var.cx*cmp_mem;

% update mineral systems oxide compositions for solid assemblage
PHS_oxdfit = zeros(np,cal.nmsy+1,cal.noxd);
PHS_oxdfit(:,1,:) = MLT_oxdfit;
for j = 1:cal.nmsy
    PHS_oxdfit(:,j+1,:) = SOL_memfit(:,cal.msy_mem(j,:)==1)*cal.mem_oxd(cal.msy_mem(j,:)==1,:)./sum(SOL_memfit(:,cal.msy_mem(j,:)==1)+1e-16,2);
end
PHS_oxdfit = PHS_oxdfit.*(PHS>0);

% get fitted phase fractions
PHS_frcfit = zeros(size(PHS));
PHS_frcfit(:,1) = var.m*100;
PHS_frcfit(:,2:end) = SOL_memfit*cal.msy_mem.';%.*(1-PHS_frcfit(:,1)/100);

datafit = [MLT_oxdfit(:);SOL_oxdfit(:);SOL_memfit(:);PHS_frcfit(:)];%repmat(PHSfit(:,1),6,1)];

end