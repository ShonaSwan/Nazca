% Chebyshev-like fixed-point iterative update with Anderson acceleration

% res      : preconditioned residual of governing equation
% x        : current iterate
% g        : fixed-point updated iterate
% f        : current fixed-point update f = g - x
% x_acc    : fully accelerated iterate
% x_new    : new updated iterate with (some/no) acceleration applied
% rho.est  : estimated spectral radius at each grid point
% rho.mean : estimated spectral radius as global mean
% FHST     : history of previous fixed-point updates
% itpar... : iterative parameter structure 
%      .cheb.alpha : damping for first Chebychev-like coefficient
%      .cheb.beta  : damping for second Chebychev-like coefficient
%      .cheb.gamma : damping for third Chebychev-like coefficient 
%      .anda.m     : depth for Anderson acceleration
%      .anda.mix   : mixing coefficient for Anderson acceleration
%      .anda.reg   : regularisation coefficient for Anderson acceleration
% count    : iter*step count to track history of run

function [x,XHST,FHST,rho] = iterate(x,res,rho,XHST,FHST,itpar,count)

% allocate arrays of correct shape
alpha = 0.*x;
beta  = 0.*x;
gamma = 0.*x;
x_new = 0.*x;
x_acc = 0.*x;
f     = 0.*x;
g     = 0.*x;


% 1) Chebyshev-like fixed-point iterative update

% Update per-DOF spectral radius rho estimates from ratio of consecutive updates
ratio   = abs(FHST(:,end)) ./ abs(FHST(:,end-1) + 1e-12);  % form ratio of two most recent updates
rho_new = min(0.99, max(0.01, ratio));                     % clamp to [0.01, 0.99]

% Moving average for stability
rho.est  = 0.7*rho.est  + 0.3*rho_new;           % moving average
rho.mean = 0.9*rho.mean + 0.1*mean(rho.est);     % moving average
rho.est  = max(rho.est, 0.5*rho.mean);           % avoid crazy outliers
rho.est  = max(rho.est, 0.1);                    % global lower bound

% Chebyshev-like coefficients
alpha(:) =  4./(2 + rho.est).^2              .*itpar.cheb.alpha;
beta (:) =  2.*(2 - rho.est) ./ (2 + rho.est).*itpar.cheb.beta;
gamma(:) = -1./(2 + rho.est).^2              .*itpar.cheb.gamma;

% New fixed-point update and updated iterate
f(:) = -alpha(:).*res(:) + beta(:).*FHST(:,end) + gamma(:).*FHST(:,end-1);
g    = x + f;


% 2) Anderson acceleration (update-based)

% Shift histories and store current g, f
FHST = [FHST(:,2:end) f(:)];   % previous solution updates
XHST = [XHST(:,2:end) g(:)];   % previous solution updates

if count>itpar.anda.m && itpar.anda.mix>eps  % only if enough history and mix>0

    % Take differences of updates
    DF  = FHST(:,2:end) - FHST(:,1:end-1);   % n×m
    DX  = XHST(:,2:end) - XHST(:,1:end-1);   % n×m
    reg = itpar.anda.reg.*rms(DF.'*DF,'all');

    % Solve min_gamma || f - DF * delta + reg*I ||_2  (global regularised least squares)
    delta = (DF.'*DF + reg*eye(itpar.anda.m)) \ (DF.'*f(:));

    % Standard Anderson Type-I update for fixed-point:
    % x_acc(:) = g(:) - FHST(:,1:end-1) * delta;
    x_acc(:) = g(:) - DX * delta;

    % Damped Anderson step
    x_new(:) = g(:) + itpar.anda.mix * (x_acc(:) - g(:));

else
    % No acceleration applied
    x_new = g;
end

% 3) Update the current iterate with the new fixed-point/accelerated value
x = x_new;

end