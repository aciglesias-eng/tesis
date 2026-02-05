function [nsr,Fi,nfield] = solarReceiver(TIT,Edni)

alfa =0.95;       % receiver solar absorptance,
epsilon =0.85;    % emittance
fi = 5.67e-8;     % Boltzman W/m2 K4
Fview = 1;        % Radiative view factor
Tr = TIT + 150;   % surface temperature of the solar receiver 
fcon= 1;          % Convective heat loss facto8r
hconv= 10;        % W/m2 K
Tamb= 298.15;     % K
C=900;            % Concentration ratio
nfield= 0.6;      % heliostat field efficiency
% Edni=800;         % W/m2
Ts=5800;
T0=298.15;

nsr=alfa- (epsilon * fi *Fview * Tr^4  + fcon * hconv*(Tr-Tamb) )/(nfield*  Edni* C);
Fi= 1- (3/4)*(T0/Ts)*( 1-cos(0.005) )^(1/4) + (1/3)*(T0/Ts)^4;
if Edni ==0 
    nsr=0;
elseif nsr<0
    nsr=0;
end


