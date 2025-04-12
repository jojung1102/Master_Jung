function out = Bending_system_2D

out{1} = @init;
out{2} = @fuN_e_valval;
% out{3} = @jacobian;
% out{4} = @jacobianp;
% out{5} = @hessians;
% out{6} = @hessiansp;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function [tspan, X0, options] = init
handles = feval(Bending_system_2D);
X0=[0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 100];

% --------------------------------------------------------------------------

% Dynamical functions
function dXdt = fuN_e_valval(t, kmrgd, par_C0, par_zeta, par_kappa)

    A_I = kmrgd(1);
    A_cont = kmrgd(2);
    
    % Building blocks
    C = @(A) sqrt((1 - ((2.* A_cont)./A - 1).^2))./(2*sqrt(A_cont));
    C_dA = @(A) (2 .* A_cont-A) ./ ...
        (2 * A.^2 .* sqrt(A - A_cont));
    C_dA_cont = @(A) - 1 ./ (2 * A .* sqrt(A - A_cont));
    
    dA_Idt = - 2 .* par_kappa .* ((C(A_I) - par_C0).^2 - (C((1 - A_I)) - par_C0).^2 + ...
        2 .* A_I .* (C(A_I) - par_C0) .* C_dA(A_I) - 2 .* (1 - A_I) .* (C((1 - A_I)) - par_C0) .* C_dA((1 - A_I)));
    
    dA_contdt = - (4 * par_kappa .* (A_I .* (C(A_I) - par_C0) .* C_dA_cont(A_I) + ...
        (1 - A_I) .* (C((1 - A_I)) - par_C0) .* C_dA_cont((1 - A_I))) +  par_zeta*sqrt(pi)./sqrt(A_cont));
    
    dXdt = [dA_Idt; dA_contdt];


% --------------------------------------------------------------------------
function jac = jacobian(t, kmrgd, par_C0, par_zeta, par_kappa)
% --------------------------------------------------------------------------
function jacp = jacobianp(t, kmrgd, par_C0, par_zeta, par_kappa)
% --------------------------------------------------------------------------
function hess = hessians(t, kmrgd, par_C0, par_zeta, par_kappa)
% --------------------------------------------------------------------------
function hessp = hessiansp(t, kmrgd, par_C0, par_zeta, par_kappa)
%---------------------------------------------------------------------------
function tens3  = der3(t, kmrgd, par_C0, par_zeta, par_kappa)
%---------------------------------------------------------------------------
function tens4  = der4(t, kmrgd, par_C0, par_zeta, par_kappa)
%---------------------------------------------------------------------------
function tens5  = der5(t, kmrgd, par_C0, par_zeta, par_kappa)