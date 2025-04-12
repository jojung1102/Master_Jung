function out = Bending_system_4D

out{1} = @init;
out{2} = @fuN_eval;
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
function [tspan, X0,options] = init
handles = feval(Bending_system_4D);
X0=[0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 100];

% --------------------------------------------------------------------------

% Dynamical functions
function dXdt = fuN_eval(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon)

    A_I = kmrgd(1);
    A_cont = kmrgd(2);
    n_fI = kmrgd(3);
    n_dI = kmrgd(4);

    gamma = 1 ./ (1 + exp(-par_epsilon));
    
    C = @(A) sqrt((1 - ((2.* A_cont)./A - 1).^2))./(2*sqrt(A_cont));
    C_dA = @(A) (2 .* A_cont-A) ./ ...
        (2 * A.^2 .* sqrt(A - A_cont));
    C_dA_cont = @(A) - 1 ./ (2 * A .* sqrt(A - A_cont));
    
    n_I = n_dI + n_fI;
    n_II = par_n_f + par_n_d - n_I;
    C_I = par_c * n_I./A_I;
    C_II = par_c * n_II./(1 - A_I);
    
    entropic_contribution = log((A_I .* ((1 - A_I) - n_II)) ./ ((A_I - n_I) .* (1 - A_I)));
    
    % Change in variables; the free energy is rescaled by a/A
    dA_Idt = entropic_contribution - 2 .* par_kappa .* ((C(A_I) - par_C0 - C_I).^2 - (C((1 - A_I)) - par_C0 - C_II).^2 + ...
        2 .* A_I .* (C(A_I) - par_C0 - C_I) .* (C_dA(A_I)+C_I/A_I) - 2 .* (1 - A_I) .* (C((1 - A_I)) - par_C0 - C_II) .* (C_dA(1 - A_I)+C_II/(1 - A_I)));
    
    dA_contdt = - (4 * par_kappa .* (A_I .* (C(A_I) - par_C0 - C_I) .* C_dA_cont(A_I) + ...
        (1 - A_I) .* (C((1 - A_I)) - par_C0 - C_II) .* C_dA_cont((1 - A_I))) +  par_zeta*sqrt(pi)./sqrt(A_cont));
    
    dn_fIdt = -log((gamma .* n_fI .* (1 - A_I-n_II)) ./ ((par_n_f - n_fI) .* (A_I-n_I))) + 4 * par_kappa * par_c * ((C(A_I) - par_C0 - C_I) - (C(1 - A_I) - par_C0 - C_II));
    dn_dIdt = -log((n_dI .* (1 - A_I-n_II)) ./ ((par_n_d - n_dI) .* (A_I-n_I))) + 4 * par_kappa * par_c * ((C(A_I) - par_C0 - C_I) - (C(1 - A_I) - par_C0 - C_II));

    dXdt = [dA_Idt; dA_contdt; dn_fIdt; dn_dIdt]; 

% --------------------------------------------------------------------------
function jac = jacobian(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon)
% --------------------------------------------------------------------------
function jacp = jacobianp(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon)
% --------------------------------------------------------------------------
function hess = hessians(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon)
% --------------------------------------------------------------------------
function hessp = hessiansp(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon)
%---------------------------------------------------------------------------
function tens3  = der3(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon)
%---------------------------------------------------------------------------
function tens4  = der4(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon)
%---------------------------------------------------------------------------
function tens5  = der5(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon)