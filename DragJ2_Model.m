function [a, e, i, RAAN, w, M, t, lifetime] = DragJ2_Model(a0, e0, i0, RAAN0, w0, M0, period, tspan, Cd, A, m, reentry_alt, perturbation_type)
%=========================================================================%
% FUNCTION: DragJ2_Model.m
% AUTHOR: Dante Sanaei (dsanaei@purdue.edu), 2022
% DESCRIPTION: Develops a dynamic model using Guass's Variational
% Equations that incorporates perturbative forces of drag and Earth's J2
% effect. 
% INPUTS:
%   > a0 = Initial semi-major axis (km)
%   > e0 = Initial eccentricity
%   > i0 = Initial inclination (rad)
%   > w0 = Initial argument of perigee (rad)
%   > RAAN0 = Initial right ascention of ascending node (rad)
%   > M0 = Initial mean anomaly (rad)
%   > period = Initial period (s)
%   > tspan = Timespan array 
%   > Cd = Coefficient of drag of the spacecraft
%   > A = Frontal area of the spacecraft (m^2)
%   > m = Mass of spacecraft (kg)
% OUPUTS:
%   > a = semi-major axis array (km)
%   > e = eccentricity array
%   > i = inclination array (rad)
%   > w = argument of perigee array (rad)
%   > RAAN = right ascention of ascending node array (rad)
%   > M = mean anomaly array (rad)
%   > t = period array (s)
%   > lifetime = lifetime of orbit (days). nan if orbit never decays 
%=========================================================================%
    %% Perform numerical integration
    global Re

    % Init
    init = [a0 e0 i0 RAAN0 w0 M0]; 
    Re = 6378; % Radius of earth (km)

    % Run ode45 
    options = odeset('reltol',1e-8,'abstol',1e-8,'initialstep',period/10000, 'Events', @altCheck);
    [t,Y] = ode89(@integrator, tspan, init', options);

    % Outputs
    a = Y(:,1); e = Y(:,2); i = Y(:,3); RAAN = Y(:,4); w = Y(:,5); M = Y(:,6);

    % Numerical Lifetime Calculation
    lifetime = orbit_lifetime(a, e, M, t);

    %% Numerical Integration Function 
    function dGVEdt = integrator(t,F)
    % Constants
    mu = 398600; % Earths gravitational parameter
    J2 = 1082.63e-6; % Earth's J2 constant
    
    % Integration Variables:
    % semi-major axis, eccentiricty, inclination, Right Ascention of Ascending
    % Node, Argument of perigee, Mean anomaly
    a = F(1); e = F(2); i = F(3);  RAAN = F(4);  w = F(5);  M = F(6);
    
    % Extra orbital parameters
    n = sqrt(mu/a^3); % Mean motion
    f = M + (2*e-(1/4)*e^3)*sin(M) + (5/4) * e^2 * sin(2*M) + (13/12) * e^3 *sin(3*M); % True anomaly
    r = (a*(1-e^2))/(1+e*cos(f)); % radius (km)
    u = w + f; % argument of latitude
    h = sqrt(mu*a*(1-e^2)); % angular momentum
    V = sqrt( (1 + e^2 + 2*e*cos(f)) * (h^2 / (a^2 * (1-e^2)^2)) ); % velocity (km/s)
    b = (Cd*A)/(2*m); % balistic coefficient

    % Atmosphere calculations
    h_geometric = r - Re; % geometric altitude (km)
    rho = atmosphere(h_geometric); % atmosphere model
    rho = rho * (1e9); % density (kg/km^3)

    % Components of Forces in RSW frame
    if perturbation_type == perturbation.J2 || perturbation_type == perturbation.All
        R_j2 = - ((1.5*J2*mu*Re^2)/(r^4)) * (1-3* (sin(u)^2*sin(i)^2));
        S_j2 = - ((1.5*J2*mu*Re^2)/(r^4)) * sin(i)^2 * sin(2*u);
        W_j2 = - ((1.5*J2*mu*Re^2)/(r^4)) * sin(u) * sin(2*i);
    else
        R_j2 = 0; S_j2 = 0; W_j2 = 0;
    end
    
    if perturbation_type == perturbation.Drag || perturbation_type == perturbation.All
        R_drag = ((e*sin(f))/(1+e^2+2*e*cos(f))^(1/2)) * (-b*rho*V^2);
        S_drag = ((1+e*cos(f)) / (1+e^2+2*e*cos(f))^(1/2)) * (-b*rho*V^2);
        W_drag = 0;
    else
        R_drag = 0; S_drag = 0; W_drag = 0;
    end

    R = R_drag + R_j2;
    S = S_drag + S_j2;
    W = W_drag + W_j2;

    % Gauss's Form of the Variational Equations:
    dadt = ((2*e*sin(f))/(n*sqrt(1-e^2))) * R + ((2*a*sqrt(1-e^2))/(r*n)) * S;
    dedt = ((sqrt(1-e^2)*sin(f))/(a*n)) * R + ((sqrt(1-e^2))/(a^2*n*e)) * ((a^2*(1-e^2)-r^2)/r) * S;
    didt = ((r*cos(u))/(a^2*n*sqrt(1-e^2))) * W;
    dRAANdt = ((r*sin(u))/(a^2*n*sqrt(1-e^2)*sin(i))) * W;
    dwdt =  ((-sqrt(1-e^2)*cos(f))/(a*n*e)) * R + ((sqrt(1-e^2)*sin(f))/(a*n*e))*(1 + r/(a*(1-e^2))) * S - ((r*sin(u)*cot(i))/(a^2*n*sqrt(1-e^2))) * W;
    dsdt = ((((1-e^2)*cos(f))/(a*n*e)) - ((2*r)/(a^2*n))) * R - (((1-e^2)*sin(f))/(a*n*e))*(1 + r/(a*(1-e^2))) * S;
    dMdt = n + dsdt;

    % Return 
    dGVEdt = [dadt dedt didt dRAANdt dwdt dMdt]';
    end

    %% Stop ode when entering atmosphere (alt < 105km)
    function [value, isterminal, direction] = altCheck(t, Y)
        a = Y(1); e = Y(2); M = Y(6);
        f = M + (2*e-(1/4)*e^3)*sin(M) + (5/4) * e^2 * sin(2*M) + (13/12) * e^3 *sin(3*M); % True anomaly
        r = (a*(1-e^2))/(1+e*cos(f)); % Radius (km)
        h_geometric = r - Re; % geometric altitude (km)
        value      = (h_geometric <= reentry_alt);
        isterminal = 1;   
        direction  = 0;
    end

    %% Numerical lifetime of orbit
    function lifetime = orbit_lifetime(a, e, M, t)
        f = M + (2*e-(1/4)*e.^3).*sin(M) + (5/4) * e.^2 .* sin(2*M) + (13/12) * e.^3 .*sin(3*M); % True anomaly
        r = (a.*(1-e.^2))./(1+e.*cos(f)); % Radius (km)
        altitude = r - Re;  % geometric altitude (km)
        [min_alt, min_alt_index] = min(altitude);
        if min_alt <= reentry_alt
            lifetime = t(min_alt_index) / (36 * 2400); % days
        else
            lifetime = nan;
        end
    end
end