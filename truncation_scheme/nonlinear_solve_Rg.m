clear;
close all;
global vars rel mob

%% Draw optimization landscape
% vars.N = 100;
% avec = linspace(0.01,0.2,20);
% Rg = linspace(0.8,1,40);
% for j = 1:length(avec)
%     for i = 1:length(Rg)
%         vars.scale = avec(j);
%         A(i,j) = max(effective_sphere_diff(Rg(i),vars,1,1));
%     end
% end

%%
% figure()
% imagesc(avec,Rg,log10(abs(A)));



%vars = {};
mob = 1;
rel = 1;

A = [];
b = [];
Aeq = [];
beq = [];

lb = [0.01 1e-4];
ub = [1 1];

% lb = [];
% ub = [];

vars.scale = 1; 
x0 = [0.99; 0.1];

options = optimset('TolCon',1e-12,'TolFun',1e-12);

%solve optimisation problems for the different resolutions. 
Nlist = [100,700];
for i = 1:length(Nlist)
    N = Nlist(i);
    vars.N = N;
    [x,fval] = fmincon(@objfun,x0,A,b,Aeq,beq,lb,ub,@confuneq,options);
    % x = [R,L];
    bb(i,:) = x; 
    err(i) = fval;
    vars.scale =x(2);
    err_all_mob(i,:) = -effective_sphere_diff(x(1),vars,1,0);
    err_all_res(i,:) = -effective_sphere_diff(x(1),vars,0,0);
end

%%


function [c,ceq] = confuneq(x)
    global vars mob rel
    % Nonlinear inequality constraints  
    vars.scale = x(2); 
    c = -effective_sphere_diff(x(1),vars,mob,rel); 
    % Nonlinear equality constraints
    ceq = [];
end

function f = objfun(x)
    global vars rel mob  
    %if measuring the relative error instead
    vars.scale = x(2); 
    f = max(abs(effective_sphere_diff(x(1),vars,mob,rel)));
end