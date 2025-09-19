% create manufactured solution
load ocean;  
clear x z SOL W U wm um Pf Pc
TINY = 1e-16;
mumin = 1e-3;
syms U_mms(x,z) W_mms(x,z) wm_mms(x,z) um_mms(x,z) Pf_mms(x,z) Pc_mms(x,z) eta_mms(x,z)  zeta_mms(x,z) rho_mms(x,z) rhom_mms(x,z) src_mms(x,z) Ks_mms(x,z) mu_mms(x,z)

fprintf(1,'\n\n  ***  compose manufactured solution\n\n');


%Making a Gaussian for the melt fraction 
x0 = L/2;
z0 = D/2;
width = max(L,D) / 6; % 1/6 width of domain


% compose manufactured material coefficients and volume source
gp_mms(x,z) = exp(-((x - x0)^2 + (z - z0)^2) / width^2);
mu_mms(x,z) = 0.05 * gp_mms(x,z);

eta_mms(x,z)  = 1e+19 + 1e+18*(cos(4*(x)*pi/L)*sin(4*(z)*pi/L)) * exp(-3e+1 * mu_mms(x,z));
rho_mms(x,z)  = 3e+3+5e+1*(cos(4*(x)*pi/L)*sin(4*(z)*pi/L)); rhoref = 3000; %int(rho_mms(x,z),x ,0,L)/L;
rhom_mms(x,z) = 28e+2-5e+1*(cos(4*(x)*pi/L)*sin(4*(z)*pi/L)); 
src_mms(x,z)  = -4e-14*(cos(4*(x)*pi/L)*sin(4*(z)*pi/L)) * gp_mms(x,z);
zeta_mms(x,z) = eta_mms(x,z)/(mu_mms(x,z) + 1e-3); 
Ks_mms(x,z)   = mu_mms(x,z)^2 * 1e-9;
M_mms(x,z)    = mu_mms(x,z) * rhom_mms(x,z);

% compose manufactured solution variables
W_mms(x,z)  = 6e-9.*(cos(4*(x)*pi/L).*sin(4.5*(z)*pi/L)) * (1 + 3*gp_mms(x,z))/4;
U_mms(x,z)  = 5e-9.*(sin(4*(x)*pi/L).*cos(  4*(z)*pi/L)) * (1 + 3*gp_mms(x,z))/4;
DivV_mms(x,z) = (diff(W_mms,z) + diff(U_mms,x));
Pf_mms(x,z) = -1e6 .*(cos(4*(x)*pi/L).*sin(4*(z)*pi/L)) + zeta_mms(x,z) * DivV_mms(x,z);
wm_mms(x,z) = - Ks_mms(x,z) * (diff(Pf_mms, z) - ( rhom_mms(x,z) - rhoref ) * g0); %3e-10.*(sin(4*(x)*pi/L).*cos(  4*(z)*pi/L)); 
um_mms(x,z) = - Ks_mms(x,z) *  diff(Pf_mms, x); %1e-9.*(cos(4*(x)*pi/L).*sin(4.5*(z)*pi/L)); 
Pc_mms(x,z) = - DivV_mms(x,z) * zeta_mms(x,z); %5e8 .*(sin(  4*(z)*pi/L))* mu_mms(x,z);

fprintf(1,'       W    = %s \n',char(W_mms));
fprintf(1,'       U    = %s \n',char(U_mms));
fprintf(1,'       wm   = %s \n',char(wm_mms));
fprintf(1,'       um   = %s \n',char(um_mms));
fprintf(1,'       Pf   = %s \n',char(Pf_mms));
fprintf(1,'       Pc   = %s \n',char(Pc_mms));
fprintf(1,'       Mu   = %s \n',char(mu_mms));
fprintf(1,'       eta  = %s \n',char(eta_mms));
fprintf(1,'       zeta = %s \n',char(zeta_mms));
fprintf(1,'       Ks   = %s \n',char(Ks_mms));
fprintf(1,'       rho  = %s \n',char(rho_mms));
fprintf(1,'       rhom  = %s \n',char(rhom_mms));
fprintf(1,'       src  = %s \n',char(src_mms));
fprintf(1,'       . ');

% update strain rates
DivV_mms(x,z)    = (diff(W_mms,z) + diff(U_mms,x));
DivrhoV_mms(x,z) = (diff(rho_mms*W_mms,z) + diff(rho_mms*U_mms,x));
DivMv_mms(x,z)   = (diff(M_mms*wm_mms,z) + diff(M_mms*um_mms,x));
exx_mms(x,z)     = diff(U_mms,x) - DivV_mms./3;         % x-normal strain rate
ezz_mms(x,z)     = diff(W_mms,z) - DivV_mms./3;         % z-normal strain rate
exz_mms(x,z)     = 1/2*(diff(U_mms,z)+diff(W_mms,x));  % xz-shear strain rate
fprintf(1,' . ');

% update stresses
txx_mms(x,z) = eta_mms * exx_mms;                  % x-normal stress
tzz_mms(x,z) = eta_mms * ezz_mms;                  % z-normal stress
txz_mms(x,z) = eta_mms * exz_mms;                  % xz-shear stress
fprintf(1,' . ');

% manufactured solution residuals
res_W_mms = -(diff(tzz_mms(x,z),z) + diff(txz_mms(x,z),x)) + diff(Pf_mms(x,z),z) + diff(Pc_mms(x,z),z) - (rho_mms(x,z)-rhoref)*g0;
res_U_mms = -(diff(txx_mms(x,z),x) + diff(txz_mms(x,z),z)) + diff(Pf_mms(x,z),x) + diff(Pc_mms(x,z),x);
res_wm_mms =  wm_mms(x,z) + Ks_mms(x,z) * (diff(Pf_mms(x,z), z) - ( rhom_mms(x,z) - rhoref) * g0);
res_um_mms =  um_mms(x,z) + Ks_mms(x,z) *  diff(Pf_mms(x,z), x); 
res_Pf_mms =  DivrhoV_mms(x,z) + DivMv_mms(x,z) + src_mms(x,z);
res_Pc_mms =  DivV_mms(x,z)  + Pc_mms(x,z) / zeta_mms(x,z);
fprintf(1,' . ');

% plot manufactured solution
figure(15);
colormap(ocean);
subplot(2,3,1); fcontour( -W_mms*yr   ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufcat. $W$ [m/yr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,2); fcontour(  U_mms*yr   ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $U$ [m/yr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,3); fcontour( -wm_mms*yr  ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufcat. $wm$ [m/yr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,4); fcontour(  um_mms*yr  ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $um$ [m/yr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,5); fcontour(  Pf_mms/1e6 ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $Pf$ [MPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,6); fcontour(  Pc_mms/1e6 ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $Pc$ [MPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')

drawnow;
fprintf(1,' . ');

figure(16);
colormap(ocean);
subplot(2,3,1); fcontour(       rho_mms ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $\rho$ [kg/m$^3$]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,2); fcontour(log10(eta_mms),[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $\eta$ [log$_{10}$ Pas]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,3); fcontour(       src_mms ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $\dot{Volsrc} [1/s]$','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,4); fcontour(log10(zeta_mms),[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $\zeta$ [log$_{10}$ Pas]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,5); fcontour(       (Ks_mms),[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $K_s$ []','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,6); fcontour(       (mu_mms),[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $\mu$ []','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
drawnow;
fprintf(1,' . \n');


% evaluate mms source terms on appropriate coordinate grids
fprintf(1,'\n  ***  evaluate manufactured solution\n\n');
x_mms  = -h/2:h:L+h/2;
z_mms  = -h/2:h:L+h/2;
xu_mms = (x_mms(1:end-1)+x_mms(2:end))./2;
zw_mms = (z_mms(1:end-1)+z_mms(2:end))./2;

fprintf(1,'       Make it so…\n');
fprintf(1,'       . ');


% W-grid
[x,z] = meshgrid(x_mms,zw_mms);
src_W_mms = double(subs(res_W_mms)); fprintf(1,' . ');
% src_wm_mms = double(subs(res_wm_mms)); fprintf(1,' . ');
% U-grid 
[x,z] = meshgrid(xu_mms,z_mms);
src_U_mms = double(subs(res_U_mms)); fprintf(1,' . ');
% src_um_mms = double(subs(res_um_mms)); fprintf(1,' . ');
% Cell Centers 
[x,z] = meshgrid(x_mms,z_mms);
src_Pf_mms = double(subs(res_Pf_mms)); fprintf(1,' . ');
% src_Pc_mms = double(subs(res_Pc_mms)); fprintf(1,' . ');

% plot manufactured residuals and evaluated source terms
figure(17);
colormap(ocean);
subplot(2,3,1); fcontour(-res_W_mms,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar; box on; title('manufactured $W$-res','Interpreter','latex');
subplot(2,3,2); fcontour(-res_U_mms,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar; box on; title('manufactured $U$-res','Interpreter','latex');
subplot(2,3,3); fcontour(-res_Pf_mms,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar; box on; title('manufactured $P_f$-res','Interpreter','latex');
subplot(2,3,4); imagesc(x_mms,zw_mms,-src_W_mms); axis ij equal tight; colorbar; box on; title('evaluated $W$-res','Interpreter','latex');
subplot(2,3,5); imagesc(xu_mms,z_mms,-src_U_mms); axis ij equal tight; colorbar; box on; title('evaluated $U$-res','Interpreter','latex');
subplot(2,3,6); imagesc(x_mms,z_mms,src_Pf_mms); axis ij equal tight; colorbar; box on; title('evaluated $P_f$-res','Interpreter','latex');
drawnow;

% figure(18);
% colormap(ocean);
% subplot(2,3,1);fcontour(-res_um_mms,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar; box on; title('manufactured $um$-res','Interpreter','latex');
% subplot(2,3,2); fcontour(-res_Pf_mms,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar; box on; title('manufactured $Pf$-res','Interpreter','latex');
% subplot(2,3,3); fcontour(-res_Pc_mms,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar; box on; title('manufactured $Pc$-res','Interpreter','latex');
% subplot(2,3,4); imagesc(xu_mms,z_mms,-src_um_mms); axis ij equal tight; colorbar; box on; title('evaluated $um$-res','Interpreter','latex');
% subplot(2,3,5); imagesc(x_mms ,z_mms,-src_Pf_mms); axis ij equal tight; colorbar; box on; title('evaluated $Pf$-res','Interpreter','latex');
% subplot(2,3,6); imagesc(x_mms ,z_mms,-src_Pc_mms); axis ij equal tight; colorbar; box on; title('evaluated $Pc$-res','Interpreter','latex');
% drawnow;

% evaluate analytical solution on appropriate coordinate grids
%W grid 
[x,z] = meshgrid(x_mms,zw_mms);
W_mms   = double(subs(W_mms)); fprintf(1,' . ');
wm_mms  = double(subs(wm_mms)); fprintf(1,' . ');
Mz      = double(subs(M_mms)); fprintf(1,' . ');
Mz      = Mz(:,2:end-1);
rhow    = double(subs(rho_mms)); fprintf(1,' . ');
rhow    = rhow(:,2:end-1);
rhomw   = double(subs(rhom_mms)); fprintf(1,' . ');
rhomw   = rhomw(:,2:end-1);
muw_mms = double(subs(mu_mms));  
muw_mms = muw_mms(:,2:end-1);  
Ksw     = double(subs(Ks_mms)); fprintf(1,' . ');
Ksw     = Ksw(:,2:end-1);
%U grid 
[x,z] = meshgrid(xu_mms,z_mms);
U_mms   = double(subs(U_mms)); fprintf(1,' . ');
um_mms  = double(subs(um_mms)); fprintf(1,' . ');
Mx      = double(subs(M_mms)); fprintf(1,' . ');
Mx      = Mx(2:end-1,:);
rhou    = double(subs(rho_mms)); fprintf(1,' . ');
rhou    = rhou(2:end-1,:);
muu_mms = double(subs(mu_mms)); 
muu_mms = muu_mms(2:end-1,:);
Ksu     = double(subs(Ks_mms)); fprintf(1,' . ');
Ksu     = Ksu(2:end-1,:);
%Center
[x,z] = meshgrid(x_mms,z_mms);
Pf_mms  = double(subs(Pf_mms)); fprintf(1,' . ');
Pc_mms  = double(subs(Pc_mms)); fprintf(1,' . ');
eta     = double(subs(eta_mms)); fprintf(1,' . ');
eta     = eta(2:end-1,2:end-1);
zeta    = double(subs(zeta_mms)); fprintf(1,' . ');
zeta    = zeta(2:end-1,2:end-1);
Ks      = double(subs(Ks_mms)); fprintf(1,' . ');
Ks      = Ks(2:end-1,2:end-1);
rho     = double(subs(rho_mms)); fprintf(1,' . ');
rho     = rho(2:end-1,2:end-1);
VolSrc  = double(subs(src_mms)); fprintf(1,' . ');
VolSrc  = VolSrc(2:end-1,2:end-1);
%Corner 
[x,z]   = meshgrid(xu_mms,zw_mms);
etaco   = double(subs(eta_mms)); fprintf(1,' . ');

Drho   = rhow -mean(rhow,2);
Drhow  = rhow -mean(rhow,2);
Drhomw = rhomw-mean(rhow,2);

% get mapping arrays
Nx = length(x_mms)-2;
Nz = length(z_mms)-2;

NP = (Nz+2) * (Nx+2);
NW = (Nz+1) * (Nx+2);
NU = (Nz+2) * (Nx+1);
NV = NW+NU;
MapP = reshape(1:NP,Nz+2,Nx+2);
MapW = reshape(1:NW,Nz+1,Nx+2);
MapU = reshape(1:NU,Nz+2,Nx+1) + NW;

% set time stepping parameters
a1 = 1; a2 = 1; a3 = 0;
b1 = 1; b2 = 0; b3 = 0;

% set boundary conditions to free slip
sds = -1;
top_cnv = -1;
bot_cnv = -1;
BCA = {'',''};

% set ghosted index arrays
icx = [Nx,1:Nx,1];
icz = [Nz,1:Nz,1];
ifx = [Nx,1:Nx+1,2];
ifz = [Nz,1:Nz+1,2];

% Boundary condition options:
Wtop  =  0;   Wbot  = -1;  Wleft  = -1;  Wright  = -1;
Utop  =  -1;  Ubot  = -1;  Uleft  =  0;  Uright  =  0;
qDxtop = -1;
Pall = -1;

twophs  = ones(size(mu_mms));
twophsw = ones(size(muw_mms));
twophsu = ones(size(muu_mms));

iter = 1;

fprintf(1,' . \n');
   
% call fluid mechanics solver
FMtime = 0;
fluidmech;
