
function [ O , V, Isyn ] = IntFire_integrator( T_conc_KC,varargin )
% Performs the integration to calculate the [O], Isyn and V values for all
% the KCs over given time range.

% Constants
alpha = 0.94;
beta = 0.18;
gsyn = 0.05;
Esyn = 0;
Cm = 1;
gl = 0.089; % This value calculated from the Drosophila model
El = -65;
Vthresh = -51; % this is a test value. not sure if it will work

if nargin > 1
    for jj = 1:2:length(varargin)-1
        if isa(varargin{jj+1},'double');
            eval([varargin{jj}, '=', num2str(varargin{jj+1})]);
        else
            error([varargin{jj} 'must be of type double']);
        end
    end
end

% Initial Conditions
O = zeros(length(T_conc_KC(:,1)),1);
Isyn = zeros(length(T_conc_KC(:,1)),1);
V = El*ones(length(T_conc_KC(:,1)),1);

% Loop over all time points (1ms time interval)
for tpt = 1:length(T_conc_KC(1,:))-1
    tspan = [tpt, tpt+1];

    % vectorized integration to calculate O
    soln = ode2(@(t,o) alpha.*(1-o).*T_conc_KC(:,t) - beta.*o, tspan, O(:,tpt));
    O(:,tpt+1) = soln(2,:);
    for nid = 1:length(T_conc_KC(:,1))
        if O(nid,tpt+1) > 1
            O(nid,tpt+1) = 1;
        end  
    end 
    
    % Calculate Isyn
    Isyn(:,tpt+1) = gsyn.*O(:,tpt).*(V(:,tpt) - Esyn);

    % vectorized integration to calculate V
    vsoln = ode2(@(t,v) (-gl*(v - El) - Isyn(:,t)), tspan, V(:,tpt));
    V(:,tpt+1) = vsoln(2,:);
    for nid = 1:length(T_conc_KC(:,1))
        if V(nid,tpt) == 0
            V(nid,tpt+1) = El;
        elseif V(nid,tpt+1) > Vthresh
            V(nid,tpt+1) = 0;
        end        
    end    
end 

end

