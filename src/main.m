% initialise model run
init;

%%%%  Physical Time Stepping Loop %%%%
iterglobal = 1;
while time <= tend && step <= Nt
    
    % Time step info
    timing;

    % Store previous solution
    store;

    %%%% Reset residuals and iteration count
    resnorm  = 1;
    resnorm0 = resnorm;
    iter     = 1;
    first    = 0;

    %%%%% Non-Linear Iteration Loop %%%%
    while resnorm/resnorm0 >= rtol/(1 + frst*100) && resnorm >= atol/(1 + frst*10) && iter <= maxit*(1 + frst) || iter <=3
        
        %%%% solve thermo-chemical equations
        thermochem;

        %%%% solve fluid-mechanics equations
        fluidmech;

        %%%% update non-linear parameters and auxiliary variables
        update;

        %%%% report convergence
        report;

        %%%% Increment iteration count 
        iter = iter+1; 
        iterglobal = iterglobal+1;
   
    %%%% The end of the Non-Linear Iteration Loop %%%%
    end

    %%%% Upsate phase equlibrium
    phseql; 

    % record model history
    history;

    % Extracted bounday conditions
    twophs  = double(deltac(icz,icx)>h/4);
    twophsw = double(muw(:  ,icx)>1e-9);
    twophsu = double(muu(icz,:  )>1e-9);
    
    % Cm  = cm.*m.*rho;
    % Cx  = cx.*x.*rho;
    % C   = Cm + CCx;
    % M   = sum(Cm,3);
    % X   = sum(Cx,3);
    % RHO = M+X;
    % Sm  = M.*sm;
    % Sx  = X.*sx;
    % S   = Sm + Sx;

    % print model diagnostics
    diagnose;
    
    if tracer_sw == 1
    u_tracer = interp2(XXu, ZZu, U, tracer_X, tracer_Z, 'linear');  
    w_tracer = interp2(XXw, ZZw, W, tracer_X, tracer_Z, 'linear');
    tracer_X = tracer_X + u_tracer * dt;
    tracer_Z = tracer_Z + w_tracer * dt;
    end

    % plot model results
    if ~mod(step,nop); output; end
    
    % increment time/step
    time = time+dt;
    step = step+1;
    frst = 0;
 
%%%% The end of the Physical Time Stepping Loop %%%% 
end

% save final state of model
output;

diary off
