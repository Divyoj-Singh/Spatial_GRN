function out = intercations
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t, kmrgd, alpha, alpha0, beta, n)

dydt=[-x(1) + alpha/(1+kmrgd(5)^n)+alpha0; -kmrgd(2) + alpha/(1+kmrgd(6)^n)+alpha0; -kmrgd(3) + alpha/(1+kmrgd(4)^n)+alpha0; -beta*(kmrgd(4)-kmrgd(1)); -beta*(kmrgd(5)-kmrgd(2)); -beta*(kmrgd(6)-kmrgd(3))];


% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(interactions);
y0=[0,0,0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,kmgr,kgr,gmgr,ggr,gr0mz,z0mgr,ngrmz,nzmgr,lamdagrmz,lamdazmgr)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,kmgr,kgr,gmgr,ggr,gr0mz,z0mgr,ngrmz,nzmgr,lamdagrmz,lamdazmgr)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,kmgr,kgr,gmgr,ggr,gr0mz,z0mgr,ngrmz,nzmgr,lamdagrmz,lamdazmgr)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,kmgr,kgr,gmgr,ggr,gr0mz,z0mgr,ngrmz,nzmgr,lamdagrmz,lamdazmgr)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,kmgr,kgr,gmgr,ggr,gr0mz,z0mgr,ngrmz,nzmgr,lamdagrmz,lamdazmgr)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,kmgr,kgr,gmgr,ggr,gr0mz,z0mgr,ngrmz,nzmgr,lamdagrmz,lamdazmgr)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,kmgr,kgr,gmgr,ggr,gr0mz,z0mgr,ngrmz,nzmgr,lamdagrmz,lamdazmgr)
