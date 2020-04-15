function LE_est = ET(Ta, qa, gs, ga, Rn, G, P, gg, gr)
% Estimate ET using the alternative to the Penman-Monteith equation derived
% in McColl (2020). The equation uses exactly the same inputs as the
% Penman-Monteith equation.
% 
% INPUTS
% Ta: near-surface air temperature [K]
% qa: near-surface specific humidity [-]
% gs: surface conductance [m/s]
% ga: aerodynamic conductance [m/s]
% Rn: net radiation [W m^2]
% G: ground heat flux [W m^2]
% P: near-surface air pressure [Pa]
% gg: storage conductance [m/s]
% gr: radiative conductance [m/s]
%
% OUTPUTS:
% LE_est: estimated latent heat flux [W m^2]
%
% EXAMPLES:
% 1. Given fixed observations of net radiation and ground heat flux,  
% estimate ET using the "radiatively uncoupled" equation (analogous to the
% "radiatively uncoupled" Penman-Monteith equation, using the terminology
% of Raupach (2001)) by setting gg = gr = 0:
% [LE_est, Ts_est] = ET(Ta, qa, gs, ga, Rn, G, P, 0, 0)
%
% 2. To include the surface temperature-dependence of outgoing longwave
% radiation and ground heat flux, estimate ET using the "radiatively
% coupled" equation (analogous to the "radiatively coupled" Penman-Monteith
% equation, using the terminology of Raupach (2001)) by choosing positive
% values of gg and gr. See Appendix B of McColl (2020) for more details:
% [LE_est, Ts_est] = ET(Ta, qa, gs, ga, Rn, G, P, gg, gr)
%
% NOTE:
% If gr > 0, Rn should be interpreted as Rn*, as defined in Appendix B
% of McColl (2020).
% If gg > 0, G should be interpreted as G*, as defined in Appendix B of
% McColl (2020.
%
% Please cite:
% McColl, K.A. (2020). Practical and theoretical benefits of an alternative
% to the Penman-Monteith evapotranspiration equation. Water Resources
% Research, 56. https://doi.org/10.1029/2020WR027106
%
% Other references:
% Raupach, M.R. (2001). Combination theory and equilibrium evaporation. 
% Quarterly Journal of the Royal Meteorological Society 127, 1149?1181.


% Parameters
BOLTZMAN = 1.380658e-23;            % BOLTZMAN,  Boltzman constant (J/K)
AVOGADRO = .602214199e24;           % AVOGADRO, Avogadro constant (1/mol)
MD       = 28.9644e-3;              % MD, molar mass dry air (kg/mol)
MV  = 18.0153e-3;                   % MD, molae mass water vapor (kg/mol)
r_v = (AVOGADRO)*(BOLTZMAN) / (MV); % r_v, gas constant for water vapor (J/(kg-K))
r_d = (AVOGADRO)*(BOLTZMAN) / (MD); % r_d, gas constant for dry air (J/(kg-K))
cp  = 7./2*(r_d);                   % cp, specific heat of air (J/(kg-K))
lambda  = 2.5008e6;                 % latent heat of vaporization (J/kg)
rho = 1.2;                          % air density [kg m^-3]

% Estimate saturation specific humidity [-]
qstarTa = nan(size(Ta));
for i = 1:length(Ta)
    qstarTa(i) = calc_q_sat(Ta(i), P(i));
end

% Estimate latent heat flux using equation (B5) of McColl (2020)
LE_est = (cp.*(ga+gg+gr).*rho.*r_v.*(Ta.^2)./lambda).*lambertw(0,...
    (((gs.*ga)./(gs+ga)).*(lambda^2).*qstarTa./(cp.*r_v.*(Ta.^2).*(ga+gg+gr))).*...
    exp((lambda./(cp.*rho.*r_v.*(Ta.^2).*(ga+gg+gr))).*...
    ((Rn-G)+((gs.*ga)./(gs+ga)).*lambda.*rho.*qa))) - ...
    ((gs.*ga)./(gs+ga)).*rho.*lambda.*qa;

end

