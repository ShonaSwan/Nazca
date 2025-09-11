clear; close all;

% load default parameters
run('../usr/par_default')

 % test decreasing grid step
NN = 50.*[1,2,4]; 

for nn = NN

    % set run parameters
    runID    =  'bnchm_VP';          % run identifier
    nop      =  1;                   % output frame plotted/saved every 'nop' time steps
    bnchm    =  1;                   % set flag for mms benchmark in fluidmech

    % set model domain parameters
    N        =  nn;                  % number of grid points in z-direction (incl. 2 ghosts)
    h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]

    % update inner indeces
    inz = 2:N+1;
    inx = 2:N+1;

    % create output directory
    if ~isfolder([outdir,'/',runID])
        mkdir([outdir,'/',runID]);
    end

    % run code
    run('../src/mms')

    figure(17); clf;
    colormap(ocean);
    subplot(2,3,1); imagesc(x_mms,zw_mms,-W*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $W$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,2); imagesc(xu_mms,z_mms, U*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $U$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,3); imagesc(x_mms,zw_mms,-wm*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $wm$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,4); imagesc(x_mms,zw_mms,-(W-W_mms)*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $W$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,5); imagesc(xu_mms,z_mms, (U-U_mms)*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $U$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,6); imagesc(x_mms,zw_mms,-(wm-wm_mms)*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $wm$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    drawnow;

    figure(18); clf;
    colormap(ocean);
    subplot(2,3,1); imagesc(xu_mms,z_mms, um*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $um$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,2); imagesc(x_mms ,z_mms, Pf/1e3); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $Pf$ [kPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,3); imagesc(x_mms ,z_mms, Pc/1e3); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $Pc$ [kPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,4); imagesc(xu_mms,z_mms, (um-um_mms)*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $um$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,5); imagesc(x_mms ,z_mms, (Pf-Pf_mms)/1e3); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $Pf$ [kPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,6); imagesc(x_mms ,z_mms, (Pc-Pc_mms)/1e3); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $Pc$ [kPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    drawnow;

    % get solution error
    EW = norm(W-W_mms)./norm(W_mms);
    EU = norm(U-U_mms)./norm(U_mms);
    Ewm = norm(wm-wm_mms)./norm(wm_mms);
    Eum = norm(um-um_mms)./norm(um_mms);
    EPf = norm(Pf-Pf_mms)./norm(Pf_mms);
    EPc = norm(Pc-Pc_mms)./norm(Pc_mms);

    clist = [colororder;[0 0 0]];

    % plot error convergence
    fh19 = figure(19);
    p1 = loglog(h,EW, 's','Color',clist(1,:),'MarkerSize',10,'LineWidth',2); axis xy tight; hold on; box on;
    p2 = loglog(h,EU, 'o','Color',clist(2,:),'MarkerSize',10,'LineWidth',2); axis xy tight; hold on; box on;
    p3 = loglog(h,Ewm,'v','Color',clist(3,:),'MarkerSize',10,'LineWidth',2); axis xy tight; hold on; box on;
    p4 = loglog(h,Eum,'+','Color',clist(4,:),'MarkerSize',10,'LineWidth',2); axis xy tight; hold on; box on;
    p5 = loglog(h,EPf,'*','Color',clist(5,:),'MarkerSize',10,'LineWidth',2); axis xy tight; hold on; box on;
    p6 = loglog(h,EPc,'p','Color',clist(6,:),'MarkerSize',10,'LineWidth',2); axis xy tight; hold on; box on;
    set(gca,'LineWidth',1.5,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('grid step [m]','Interpreter','latex','FontSize',15)
    ylabel('rel. numerical error [1]','Interpreter','latex','FontSize',15)
    title('Numerical convergence in space','Interpreter','latex','FontSize',18)

    if nn == NN(1)
        p7 = loglog(D./NN,mean([EW,EU,Ewm,Eum,EPf,EPc]).*(NN(1)./NN).^2,'k-','LineWidth',2);  % plot linear trend for comparison
    end
    if nn == NN(end)
        legend([p1,p2,p3,p4,p5,p6,p7],{'error W','error U','error wm','error um','error Pf','error Pc','quadratic'},'Interpreter','latex','box','on','location','southeast')
    end

    % plot error convergence
    fh20 = figure(20);
    DOFS = (NN+2).*(NN+2) + 2.*(NN+1).*(NN+2);
    dofs = (nn+2).*(nn+2) + 2.*(nn+1).*(nn+2);
    p8 = loglog(dofs,FMtime,'h','Color',clist(2,:),'MarkerSize',10,'LineWidth',2); axis xy tight; hold on; box on;
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('\# dofs [1]','Interpreter','latex','FontSize',15)
    ylabel('time to solution [s]','Interpreter','latex','FontSize',15)
    title('Scaling of direct solver','Interpreter','latex','FontSize',18)

    if nn == NN(1)
        p9 = loglog(DOFS,FMtime*(DOFS./DOFS(1)).^1,'k-','LineWidth',2);  % plot linear trend for comparison
        p10 = loglog(DOFS,FMtime*(DOFS./DOFS(1)).^2,'k--','LineWidth',2);  % plot linear trend for comparison
    end
    if nn == NN(end)
        legend([p8,p9,p10],{'time to solution','linear','quadratic'},'Interpreter','latex','box','on','location','southeast')
    end

end

name = [outdir,'/',runID,'/',runID,'_bnchm'];
print(fh19,name,'-dpng','-r300','-vector');

name = [outdir,'/',runID,'/',runID,'_sclng'];
print(fh20,name,'-dpng','-r300','-vector');
