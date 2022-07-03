function [a e i RAAN w M] = analytical_model(a0, e0, i0, RAAN0, w0, M0, tspan)
    %% Analytically solve LPE's

    % Constants
    mu = 398600;
    Re = 6378; 
    J2 = 1082.63e-6;
    
    % Analytically solved Lagrange's Planetary Equations for J2 Perturbation:
    n = sqrt(mu/a0^3);
    a = a0 * ones(1, length(tspan));
    e = e0 * ones(1, length(tspan));
    i = i0 * ones(1, length(tspan));
    RAAN = RAAN0 + (-((3*n*J2) / (2*(1-e0^2)^2)) * (Re/a0)^2 * cos(i0)) .* tspan;
    w =  w0 + (-((3*n*J2) / (2*(1-e0^2)^2))* (Re/a0)^2 * (5/2 * sin(i0)^2 - 2)) .* tspan;
    M = M0 + n * tspan + (-((3*n*J2)/(2*sqrt((1-e0^2)^3))) * (Re/a0)^2 * (3/2 * sin(i0)^2 - 1)) .* tspan;
end