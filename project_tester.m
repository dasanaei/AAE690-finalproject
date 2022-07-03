%% Definitions
%Conversion Factors
hours = 3600; %Hours to seconds
days = 24*hours; %Days to seconds
deg = pi/180; %Degrees to radians

%Constants
mu = 398600; %Gravitational parameter (km^3/s^2)
RE = 6378; %Earth’s radius (km)

%Initial orbital parameters:
rp0 = RE + 215; %perigee radius (km)
ra0 = RE + 939; %apogee radius (km)
RA0 = 340*deg; %Right ascencion of the node (radians)
i0 = 65.1*deg; %Inclination (radians)
w0 = 58*deg; %Argument of perigee (radians)
TA0 = 332*deg; %True anomaly (radians)

e0 = (ra0 - rp0)/(ra0 + rp0); %eccentricity
M0 = TA0 - 2*e0*sin(TA0) + (0.75*e0^2 + (1/8)*e0^4)*sin(2*TA0) - (1/3)*(e0^3) * sin(3*TA0) + (5/32)*(e0^4)*sin(4*TA0); %mean anomaly (rad)
h0 = sqrt(rp0*mu*(1 + e0)); %angular momentrum (km^2/s)
a0 = (rp0 + ra0)/2; %Semimajor axis (km)
T0 = 2*pi/sqrt(mu)*a0^1.5; %Period (s)

CD = 2.2; 
A = (pi/4*(1^2)) * (1e-6); 
m = 100;

%% Solution
init = [a0 e0 w0 M0 i0 RA0]; %initial orbital elements

t0 = 0;
tf = 150.0*days;
nout = 20000; %Number of solution points to output for plotting purposes
tspan = linspace(t0, tf, nout);
[a e i RA w M t lifetime] = DragJ2_Model(a0, e0, i0, RA0, w0, M0, T0, tspan, CD, A, m, perturbation.Drag);
lifetime

nu = M + (2*e-(1/4)*e.^3).*sin(M) + 5/4 * e.^2 .* sin(2*M) + 13/12 * e.^3 .* sin(3*M); % Fourier series expansion to convert MA to TA (up to third term)
r_total = [];
v_total = [];
for j = 1:length(a)
    [r, v] = keplerian2ijk(a(j), e(j), rad2deg(i(j)), wrapTo360(rad2deg(RA(j))), wrapTo360(rad2deg(w(j))), wrapTo360(rad2deg(nu(j))));
    r_total = [r_total, r];
end
figure(4)
comet3(r_total(1,:), r_total(2,:), r_total(3,:))
title("J2 Perturbed Orbit of ISS over 24 Hours (Epoch = 11 March 2022 21:22:22)")
xlabel("x position (km)")
ylabel("y position (km)")
zlabel("z position (km)")
grid on;
%% Plotting

h = sqrt(mu*a.*(1-e.^2)); %angular momentum (km^2/s)
f = M + (2.*e-(1/4).*e.^3).*sin(M) + (5/4) .* e.^2 .* sin(2.*M) + (13/12) .* e.^3 .*sin(3.*M); 
p = (1-e.^2);
a = (h.^2) ./ (mu.*(1-e.^2));
r = (a.*p)./(1+e.*cos(f)); %radius
u = w + f; %argument of latitude

xk_p = r .* cos(u);
yk_p = r .* sin(u);

xk = xk_p .* cos(RA) - yk_p .* cos(i) .* sin(RA);
yk = xk_p .* sin(RA) + yk_p .* cos(i) .* cos(RA);
zk = yk_p .* sin(i);

altitude = r-RE;

%
figure(1)
subplot(3,2,1)
plot(t/days,RA/deg)
title('RAAN (deg)')
xlabel('Time (days)')
grid on
axis tight

subplot(3,2,2)
plot(t/days,w/deg)
title('Argument of Perigee (deg)')
xlabel('Time (days)')
grid on
axis tight

subplot(3,2,3)
plot(t/days,e)
title('Eccentricity')
xlabel('Time (days)')
grid on
axis tight

subplot(3,2,4)
plot(t/days,i/deg)
title('Inclination (deg)')
xlabel('Time (days)')
grid on
axis tight

subplot(3,2,5)
plot(t/days,a)
title('Semi-Major Axis (km)')
xlabel('Time (days)')
grid on
axis tight

[max_altitude, imax] = findpeaks(altitude);
[min_altitude, imin] = findpeaks(-altitude);
maxima = [t(imax) max_altitude]; %Maximum altitudes and times
minima = [t(imin) -min_altitude]; %Minimum altitudes and times

apogee = sortrows(maxima,1); %Maxima sorted with time
perigee = sortrows(minima,1); %Minima sorted with time

figure(2)
apogee(1,2) = NaN;
plot(apogee(:,1)/days, apogee(:,2),'b','linewidth',2)
hold on
plot(perigee(:,1)/days, perigee(:,2),'r','linewidth',2)
grid on
grid minor
xlabel('Time (days)')
ylabel('Altitude (km)')
ylim([0 1000]);
legend('Apogee','Perigee');
title('Orbit Decay (Apogee/Perigee)');

figure(3)
plot(t/days, altitude);
xlabel('Time (days)');
ylabel('Altitude (km)');
title('Orbit Decay');
