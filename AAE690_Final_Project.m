%=========================================================================%
% PROGRAM: AAE690_Final_Project.m
% AUTHOR: Dante Sanaei (dsanaei@purdue.edu), 2022
% DESCRIPTION: The following analysis uses a drag+j2 perturbation model to
% study lifetimes and behavior of several popular orbit types. 
%=========================================================================%
% Init
clc; clear all; close all;
mu = 398600; % Earth's Gravitational parameter 
Re = 6378; % Earth s radius (km)
reentry_alt = 90; % Re-entry altitude (geometric) (km)

%% Part 1: Building the Model

% Satellite Data
Cd = 2.2;
m = 300; % mass (kg)
r_cross = 3; % cross sectional radius (meters)
A = (pi/4*(r_cross^2)) * (1e-6); % frontal area (km^2)

% Initial orbital elements
Rp0 = Re + 160; % perigee radius (km)
Ra0 = Re + 939; % apogee radius (km)
RAAN0 = deg2rad(339.94); % Right ascencion of the node (radians)
i0 = deg2rad(28.45); % Inclination (radians)
w0 = deg2rad(58); % Argument of perigee (radians)
f0 = deg2rad(332); % True anomaly (radians)
e0 = (Ra0-Rp0)/(Ra0+Rp0); % eccentricity
a0 = (Rp0 + Ra0)/2; % Semimajor axis (km)
M0 = f0 - 2*e0*sin(f0) + (0.75*e0^2 + (1/8)*e0^4)*sin(2*f0) - (1/3)*(e0^3) * sin(3*f0) + (5/32)*(e0^4)*sin(4*f0); % mean anomaly (rad)
Period = 2*pi/sqrt(mu)*a0^1.5; % Period (s)

% Scenario time parameters
t0 = 0;
days = 15;
tf = days * 24 * 3600; 
tspan = linspace(t0, tf, 10000);

% J2 only analysis
[a e i RAAN w M t lifetime] = DragJ2_Model(a0, e0, i0, RAAN0, w0, M0, Period, tspan, Cd, A, m, reentry_alt, perturbation.J2);
[a_av e_av i_av RAAN_av w_av M_av] = analytical_model(a0, e0, i0, RAAN0, w0, M0, tspan);

f = M + (2*e-(1/4)*e.^3).*sin(M) + 5/4 * e.^2 .* sin(2*M) + 13/12 * e.^3 .* sin(3*M); % Fourier series expansion to convert M to f (up to third term)
plotOrbit(a, e, i, RAAN, w, f, "J2 Perturbed Orbit Over 15 Days")
plotElementsandCompare(t, a, e, i, RAAN, w, f, a_av, e_av, i_av, RAAN_av, w_av, "J2 Perturbed Orbit Elements Over 15 Days")

% Drag only analysis
[a e i RAAN w M t lifetime] = DragJ2_Model(a0, e0, i0, RAAN0, w0, M0, Period, tspan, Cd, A, m, reentry_alt, perturbation.Drag);
fprintf("Part 1: Drag only lifetime -> %d \n", lifetime)
f = M + (2*e-(1/4)*e.^3).*sin(M) + 5/4 * e.^2 .* sin(2*M) + 13/12 * e.^3 .* sin(3*M); % Fourier series expansion to convert M to f (up to third term)
plotOrbit(a, e, i, RAAN, w, f, "Drag Perturbed Orbit Over 15 Days")
plotElements(t, a, e, i, RAAN, w, f, "Drag Perturbed Orbit Elements Over 15 Days")
plotAltitudes(t, a, e, i, RAAN, w, f, "Orbit Decay Over 15 Days")

% Drag + J2 analysis
[a e i RAAN w M t lifetime] = DragJ2_Model(a0, e0, i0, RAAN0, w0, M0, Period, tspan, Cd, A, m, reentry_alt, perturbation.All);
fprintf("Part 1: Drag+J2 lifetime -> %d \n", lifetime)
f = M + (2*e-(1/4)*e.^3).*sin(M) + 5/4 * e.^2 .* sin(2*M) + 13/12 * e.^3 .* sin(3*M); % Fourier series expansion to convert M to f (up to third term)
plotOrbit(a, e, i, RAAN, w, f, "Drag and J2 Perturbed Orbit Over 15 Days")
plotElements(t, a, e, i, RAAN, w, f, "Drag and J2 Perturbed Orbit Elements Over 15 Days")
plotAltitudes(t, a, e, i, RAAN, w, f, "Orbit Decay Over 15 Days")

%% Part 2a: Analysis of Sun-Synchronous Orbit 

% SSO Satellite Data
Cd = 3.3;
m = 6804; % mass (kg)
diameter = 3; % cross sectional radius (meters)
length = 5; % Length
A = diameter * length * (1e-6); % frontal area (km^2)

% Initial orbital elements 
Rp0 = Re +  650; % perigee radius (km)
Ra0 = Re + 550; % apogee radius (km)
RAAN0 = deg2rad(100); % Right ascencion of the node (radians)
e0 = (Ra0-Rp0)/(Ra0+Rp0); % eccentricity
a0 = (Rp0 + Ra0)/2; % Semimajor axis (km)

nodal_precession = 0.19910213e-6;
j2 = 1082.63e-6;
i0 = acos((-2/3) * nodal_precession * 1/j2 * ((a0*(1-e0^2))/Re)^2 * sqrt(a0^3/mu))
w0 = deg2rad(0); % Argument of perigee (radians)
f0 = deg2rad(0); % True anomaly (radians)
M0 = f0 - 2*e0*sin(f0) + (0.75*e0^2 + (1/8)*e0^4)*sin(2*f0) - (1/3)*(e0^3) * sin(3*f0) + (5/32)*(e0^4)*sin(4*f0); % mean anomaly (rad)
Period = 2*pi/sqrt(mu)*a0^1.5; % Period (s)

% Scenario time parameters
t0 = 0;
days = 50;
tf = days * 24 * 3600; 
tspan = linspace(t0, tf, 10000);

% Numerical Lifetime of Sun-Synchronous Orbit analysis
[a e i RAAN w M t lifetime] = DragJ2_Model(a0, e0, i0, RAAN0, w0, M0, Period, tspan, Cd, A, m, reentry_alt, perturbation.All);
f = M + (2*e-(1/4)*e.^3).*sin(M) + 5/4 * e.^2 .* sin(2*M) + 13/12 * e.^3 .* sin(3*M); % Fourier series expansion to convert M to f (up to third term)
plotOrbit(a, e, i, RAAN, w, f, "600 km Sun-Synchronous Orbit Over 50 Days")
plotElements(t, a, e, i, RAAN, w, f, "600 km Sun-Synchronous Orbit Elements Over 50 Days")
plotAltitudes(t, a, e, i, RAAN, w, f, "600 km Sun-Synchronous Orbit Decay")


%% Part 2b: Analysis of Lower Sun-Synchronous Orbit 

% SSO Satellite Data
Cd = 3.3;
m = 6804; % mass (kg)
diameter = 3; % cross sectional radius (meters)
length = 5; % Length
A = diameter * length * (1e-6); % frontal area (km^2)

% Initial orbital elements 
Rp0 = Re +  250; % perigee radius (km)
Ra0 = Re + 300; % apogee radius (km)
RAAN0 = deg2rad(100); % Right ascencion of the node (radians)
e0 = (Ra0-Rp0)/(Ra0+Rp0); % eccentricity
a0 = (Rp0 + Ra0)/2; % Semimajor axis (km)

nodal_precession = 0.19910213e-6;
j2 = 1082.63e-6;
i0 = acos((-2/3) * nodal_precession * 1/j2 * ((a0*(1-e0^2))/Re)^2 * sqrt(a0^3/mu))
w0 = deg2rad(0); % Argument of perigee (radians)
f0 = deg2rad(0); % True anomaly (radians)
M0 = f0 - 2*e0*sin(f0) + (0.75*e0^2 + (1/8)*e0^4)*sin(2*f0) - (1/3)*(e0^3) * sin(3*f0) + (5/32)*(e0^4)*sin(4*f0); % mean anomaly (rad)
Period = 2*pi/sqrt(mu)*a0^1.5; % Period (s)

% Scenario time parameters
t0 = 0;
days = 50;
tf = days * 24 * 3600; 
tspan = linspace(t0, tf, 10000);

% Numerical Lifetime of Sun-Synchronous Orbit analysis
[a e i RAAN w M t lifetime] = DragJ2_Model(a0, e0, i0, RAAN0, w0, M0, Period, tspan, Cd, A, m, reentry_alt, perturbation.All);
f = M + (2*e-(1/4)*e.^3).*sin(M) + 5/4 * e.^2 .* sin(2*M) + 13/12 * e.^3 .* sin(3*M); % Fourier series expansion to convert M to f (up to third term)
fprintf("Part 4: SSO lifetime -> %d \n", lifetime)
plotOrbit(a, e, i, RAAN, w, f, "275 km Sun-Synchronous Orbit Until Re-Entry")
plotElements(t, a, e, i, RAAN, w, f, "275 km Sun-Synchronous Orbit Elements Until Re-Entry")
plotAltitudes(t, a, e, i, RAAN, w, f, "275 km Sun-Synchronous Orbit Decay")

%% Part 2c: Analysis of different SSO Lifetimes
lifetimes = []; crit_alts = [];
for alt = [100 150 200 250 300 350]
    % SSO Satellite Data
    Cd = 3.3;
    m = 6804; % mass (kg)
    diameter = 3; % cross sectional radius (meters)
    length = 5; % Length
    A = diameter * length * (1e-6); % frontal area (km^2)
    
    % Initial orbital elements 
    Rp0 = Re +  alt; % perigee radius (km)
    Ra0 = Re + alt + 50; % apogee radius (km)
    RAAN0 = deg2rad(100); % Right ascencion of the node (radians)
    e0 = (Ra0-Rp0)/(Ra0+Rp0); % eccentricity
    a0 = (Rp0 + Ra0)/2; % Semimajor axis (km)
    
    nodal_precession = 0.19910213e-6;
    j2 = 1082.63e-6;
    i0 = acos((-2/3) * nodal_precession * 1/j2 * ((a0*(1-e0^2))/Re)^2 * sqrt(a0^3/mu));
    w0 = deg2rad(0); % Argument of perigee (radians)
    f0 = deg2rad(0); % True anomaly (radians)
    M0 = f0 - 2*e0*sin(f0) + (0.75*e0^2 + (1/8)*e0^4)*sin(2*f0) - (1/3)*(e0^3) * sin(3*f0) + (5/32)*(e0^4)*sin(4*f0); % mean anomaly (rad)
    Period = 2*pi/sqrt(mu)*a0^1.5; % Period (s)
    
    % Scenario time parameters
    t0 = 0;
    days = 300;
    tf = days * 24 * 3600; 
    tspan = linspace(t0, tf, 20000);
    
    % Numerical Lifetime of Sun-Synchronous Orbit analysis
    [a e i RAAN w M t lifetime] = DragJ2_Model(a0, e0, i0, RAAN0, w0, M0, Period, tspan, Cd, A, m, reentry_alt, perturbation.All);
    lifetimes = [lifetimes, lifetime];
    crit_alts = [crit_alts, alt];
end

figure
plot(crit_alts, lifetimes, '*-')
title('SSO Lifetimes From 100 km to 350 km')
xlabel('SSO Initial Altitude (km)')
ylabel('Decay Time (Days)')
grid on; 

%% Part 3a: Analysis of LEO Criticially Inclined Orbit 

% SSO Satellite Data
Cd = 3.3;
m = 6804; % mass (kg)
diameter = 3; % cross sectional radius (meters)
length = 5; % Length
A = diameter * length * (1e-6); % frontal area (km^2)

% Initial orbital elements 
Rp0 = Re +  250; % perigee radius (km)
Ra0 = Re + 300; % apogee radius (km)
RAAN0 = deg2rad(100); % Right ascencion of the node (radians)
e0 = (Ra0-Rp0)/(Ra0+Rp0); % eccentricity
a0 = (Rp0 + Ra0)/2; % Semimajor axis (km)
i0 = deg2rad(63.4);
w0 = deg2rad(270); % Argument of perigee (radians)
f0 = deg2rad(0); % True anomaly (radians)
M0 = f0 - 2*e0*sin(f0) + (0.75*e0^2 + (1/8)*e0^4)*sin(2*f0) - (1/3)*(e0^3) * sin(3*f0) + (5/32)*(e0^4)*sin(4*f0); % mean anomaly (rad)
Period = 2*pi/sqrt(mu)*a0^1.5; % Period (s)

% Scenario time parameters
t0 = 0;
days = 50;
tf = days * 24 * 3600; 
tspan = linspace(t0, tf, 10000);

% Numerical Lifetime of Sun-Synchronous Orbit analysis
[a e i RAAN w M t lifetime] = DragJ2_Model(a0, e0, i0, RAAN0, w0, M0, Period, tspan, Cd, A, m, reentry_alt, perturbation.All);
f = M + (2*e-(1/4)*e.^3).*sin(M) + 5/4 * e.^2 .* sin(2*M) + 13/12 * e.^3 .* sin(3*M); % Fourier series expansion to convert M to f (up to third term)
fprintf("Part 4: Crit Inclined lifetime -> %d \n", lifetime)
plotOrbit(a, e, i, RAAN, w, f, "250 km Critically Inclined Orbit Until Re-Entry")
plotElements(t, a, e, i, RAAN, w, f, "250 km Critically Inclined Orbit Elements Until Re-Entry")
plotAltitudes(t, a, e, i, RAAN, w, f, "250 km Critically Inclined Orbit Decay")

%% Part 3b: Analysis of different Crit Inclined Lifetimes
lifetimes = []; crit_alts = [];
for alt = [100 150 200 250 300 350]
    % SSO Satellite Data
    Cd = 3.3;
    m = 6804; % mass (kg)
    diameter = 3; % cross sectional radius (meters)
    length = 5; % Length
    A = diameter * length * (1e-6); % frontal area (km^2)
    
    % Initial orbital elements 
    Rp0 = Re +  alt; % perigee radius (km)
    Ra0 = Re + alt + 50; % apogee radius (km)
    RAAN0 = deg2rad(100); % Right ascencion of the node (radians)
    e0 = (Ra0-Rp0)/(Ra0+Rp0); % eccentricity
    a0 = (Rp0 + Ra0)/2; % Semimajor axis (km)
    i0 = deg2rad(63.4);
    w0 = deg2rad(270); % Argument of perigee (radians)
    f0 = deg2rad(0); % True anomaly (radians)
    M0 = f0 - 2*e0*sin(f0) + (0.75*e0^2 + (1/8)*e0^4)*sin(2*f0) - (1/3)*(e0^3) * sin(3*f0) + (5/32)*(e0^4)*sin(4*f0); % mean anomaly (rad)
    Period = 2*pi/sqrt(mu)*a0^1.5; % Period (s)
    
    % Scenario time parameters
    t0 = 0;
    days = 300;
    tf = days * 24 * 3600; 
    tspan = linspace(t0, tf, 10000);
    
    % Numerical Lifetime of Sun-Synchronous Orbit analysis
    [a e i RAAN w M t lifetime] = DragJ2_Model(a0, e0, i0, RAAN0, w0, M0, Period, tspan, Cd, A, m, reentry_alt, perturbation.All);
    lifetimes = [lifetimes, lifetime];
    crit_alts = [crit_alts, alt];
end

figure
plot(crit_alts, lifetimes, '*-')
title('Critically Inclined Lifetimes From 100 km to 350 km')
xlabel('Initial Altitude (km)')
ylabel('Decay Time (Days)')
grid on; 

%% Part 5: Analysis of Geostationary Transfer Orbit 

% ULA Centaur Data
Cd = 3.5;
m = 2462; % mass (kg)
diameter = 3.05; % cross sectional radius (meters)
length = 12.68; % Length
A = diameter * length * (1e-6); % frontal area (km^2)

% Initial orbital elements
Rp0 = Re +  200; % perigee radius (km)
Ra0 = Re + 35786; % apogee radius (km)
RAAN0 = deg2rad(190); % Right ascencion of the node (radians)
i0 = deg2rad(28.45); % Inclination (radians)
w0 = deg2rad(180); % Argument of perigee (radians)
f0 = deg2rad(180); % True anomaly (radians)
e0 = (Ra0-Rp0)/(Ra0+Rp0); % eccentricity
a0 = (Rp0 + Ra0)/2; % Semimajor axis (km)
M0 = f0 - 2*e0*sin(f0) + (0.75*e0^2 + (1/8)*e0^4)*sin(2*f0) - (1/3)*(e0^3) * sin(3*f0) + (5/32)*(e0^4)*sin(4*f0); % mean anomaly (rad)
Period = 2*pi/sqrt(mu)*a0^1.5; % Period (s)

% Scenario time parameters
t0 = 0;
days = 300;
tf = days * 24 * 3600; 
tspan = linspace(t0, tf, 20000);

% Numerical Lifetime of GTO analysis
[a e i RAAN w M t lifetime] = DragJ2_Model(a0, e0, i0, RAAN0, w0, M0, Period, tspan, Cd, A, m, 0, perturbation.J2);
fprintf("Part 2: GTO lifetime -> %d \n", lifetime)
f = M + (2*e-(1/4)*e.^3).*sin(M) + 5/4 * e.^2 .* sin(2*M) + 13/12 * e.^3 .* sin(3*M); % Fourier series expansion to convert M to f (up to third term)
plotOrbit(a, e, i, RAAN, w, f, "Geostationary Transfer Orbit Until Re-Entry")
plotElements(t, a, e, i, RAAN, w, f, "Geostationary Transfer Orbit Elements Until Re-Entry")
plotAltitudes(t, a, e, i, RAAN, w, f, "Geostationary Transfer Orbit Decay")

%% EXTRA: SSO Inclinations 
alt = [100:1:1000];
Rp0 = Re +  alt; % perigee radius (km)
Ra0 = Re + alt; % apogee radius (km)
RAAN0 = deg2rad(100); % Right ascencion of the node (radians)
e0 = (Ra0-Rp0)/(Ra0+Rp0); % eccentricity
a0 = (Rp0 + Ra0)/2; % Semimajor axis (km)

nodal_precession = 0.19910213e-6;
j2 = 1082.63e-6;
i0 = acosd((-2/3) * nodal_precession * 1/j2 * ((a0.*(1-e0^2))/Re).^2 .* sqrt(a0.^3./mu));

figure
plot(alt, i0)
title('Sun Synchronous Inclinations')
xlabel('Initial Orbit Altitude (km)')
ylabel('SSO Inclination (deg)')
grid on; 

%% Function to plot orbit around earth
function plotOrbit(a, e, i, RAAN, w, f, fig_title)
    r_total = [];
    for j = 1:length(a)
        [r, v] = keplerian2ijk(a(j), abs(e(j)), rad2deg(i(j)), wrapTo360(rad2deg(RAAN(j))), wrapTo360(rad2deg(w(j))), wrapTo360(rad2deg(f(j))));
        r_total = [r_total, r];
    end
    
    [fig, globe] = globe3D();
    plot3(r_total(1,:) * 1000, r_total(2,:) * 1000, r_total(3,:) * 1000)
    sgtitle(fig_title)
end

%% Function to plot orbital elements
function plotElements(t, a, e, i, RAAN, w, f, fig_title)
    days = 24*3600;
    figure
    sgtitle(fig_title)
    subplot(3,2,1)
    plot(t/days,a)
    title('Semi-Major Axis (a)')
    xlabel('Time (Days)')
    ylabel('km')
    grid on; 
    
    subplot(3,2,2)
    plot(t/days,e)
    title('Eccentricity (e)')
    xlabel('Time (Days)')
    grid on; 
    
    subplot(3,2,3)
    plot(t/days,rad2deg(i))
    title('Inclination (i)')
    xlabel('Time (Days)')
    ylabel('deg')
    grid on; 
    
    subplot(3,2,4)
    plot(t/days,rad2deg(RAAN))
    title('Right Ascension of Ascending Node (R.A.A.N.)')
    xlabel('Time (Days)')
    ylabel('deg')
    grid on; 
    
    subplot(3,2,5)
    plot(t/days,rad2deg(w))
    title('Argument of Perigee (w)')
    ylim([0, 360])
    xlabel('Time (Days)')
    ylabel('deg')
    grid on; 
end

%% Function to plot orbital decay
function plotAltitudes(t, a, e, i, RAAN, w, f, fig_title)
    days = 24*3600;
    Re = 6378; % Earth's radius (km)
    r = (a.*(1-e.^2))./(1+e.*cos(f)); % Radius (km)
    altitude = r - Re;  % geometric altitude (km)
    Rp = (a .* (1-abs(e))) - Re;
    Ra = (a .* (1+e)) - Re;
    
    figure
    sgtitle(fig_title)
    subplot(2,1,1)
    title('Orbit Altitude');
    plot(t/days, altitude);
    xlabel('Time (days)');
    ylabel('Altitude (km)');
    grid on;
    
    
    subplot(2,1,2)
    title('Apogee/Perigee Altitudes');
    plot(t/days, Ra)
    hold on
    plot(t/days, Rp)
    grid on
    xlabel('Time (days)')
    ylabel('Altitude (km)')
    legend('Apogee','Perigee');
end

%% Function to compare averaged and non-averaged orbital elements
function plotElementsandCompare(t, a, e, i, RAAN, w, f, a_av, e_av, i_av, RAAN_av, w_av, fig_title)
    days = 24*3600;
    figure
    sgtitle(fig_title)
    subplot(3,2,1)
    plot(t/days,a); hold on;
    plot(t/days,a_av); hold on;
    title('Semi-Major Axis (a)')
    xlabel('Time (Days)')
    ylabel('km')
    legend('Non-Averaged', 'Averaged')
    grid on; 
    
    subplot(3,2,2)
    plot(t/days,e); hold on;
    plot(t/days,e_av);
    title('Eccentricity (e)')
    xlabel('Time (Days)')
    legend('Non-Averaged', 'Averaged')
    grid on; 
    
    subplot(3,2,3)
    plot(t/days,rad2deg(i)); hold on
    plot(t/days,rad2deg(i_av));
    title('Inclination (i)')
    xlabel('Time (Days)')
    ylabel('deg')
    legend('Non-Averaged', 'Averaged')
    grid on; 
    
    subplot(3,2,4)
    plot(t/days,rad2deg(RAAN)); hold on
    plot(t/days,rad2deg(RAAN_av));
    title('Right Ascension of Ascending Node (R.A.A.N.)')
    xlabel('Time (Days)')
    ylabel('deg')
    legend('Non-Averaged', 'Averaged')
    grid on; 
    
    subplot(3,2,5)
    plot(t/days,rad2deg(w)); hold on
    plot(t/days,rad2deg(w_av));
    title('Argument of Perigee (w)')
    %ylim([0, 360])
    xlabel('Time (Days)')
    ylabel('deg')
    legend('Non-Averaged', 'Averaged')
    grid on; 
end










