% Inputs: 
% H = half the thickness of the tissue
% real_t = the time (in seconds) at which the data is being collected
% D = diffusion constant (in microns^2/s) of the sample
% xdata = vector of thickness points (in microns) 


function c = solvepde(D, xdata, params)

% unpack params - these are constant and will not be optimized
H = params(1); real_t = params(2); 

nondim_t = real_t/(H^2/D); % nondimensional t 

% solver parameters: don't really need to change
m = 0; 
num_tpoints = 19; 

tspan = linspace(0,1,num_tpoints);
if sum(nondim_t == tspan) == 0
    % add in data time point if doesn't already exist in array
    tspan = sort([linspace(0,1,num_tpoints) nondim_t]); 
end

% normalize xdata 
% xd = xdata/max(xdata);
% xd = xdata/H; 

% in case the data doens't have points all the way to the middle...
stepsize = xdata(2)-xdata(1);
xd = 0:stepsize/H:1; 
sol = pdepe(m, @pdefun, @icfun, @bcfun, xd, tspan);

% temp: only fit the first quarter of it
%indez = length(xdata) - floor(length(xdata)/4) + 1;
indez = length(xd)-length(xdata)+1;

% determine which index in matrix corresponds to desired time point.
t_index = tspan == nondim_t;
c = sol(t_index,indez:end,1);

end

function [c,f,s] = pdefun(x,t,u,DuDx)
c = 1;
f = DuDx;
s = 0; 
end

function u0 = icfun(x)
u0 = 0;
end

function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
pr = ur-1;
qr = 0;
pl = 0;
ql = 1;
end