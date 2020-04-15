function [ q_sat ] = calc_q_sat(T,P)
% Calculate the saturated specific humidity q_sat (kg kg^{-1}), given temperature T (K) and
% pressure P (Pa) using a standard approximation to the Clausius-Clapeyron
% relationship (valid for -30 deg C <= T <= 35 deg C)
%
% Author: Kaighin McColl
% Date: 07/11/2012

mol_ratio = 0.62198; % ratio of the molecular weights of water and dry air

T_cels = T - 273.15; % convert temperature to degrees celsius
q_sat = (mol_ratio/P)*611.2*exp((17.67*T_cels)./(T_cels + 243.5));

if T < 50
    q_sat = 0;
end


end

