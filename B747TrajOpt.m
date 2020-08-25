%% --------------------------------
%  Trajectory Optimization of B747
%  --------------------------------
%  by Hiren Saravana, March 2020
%  --------------------------------
%  Description:
%  - Nonlinear longitudinal 
%    dynamics
%  - Cruise condition at Mach 0.9
%  - Min. control effort climb in
%    altitude
%  --------------------------------

%% Main Program (Solves TPBVP Using bvp4c)

global m Jyy rho g CD S CL

altitude = 40e3;                % ft
M_inf = 0.9;
alfa = 2.4;                     % degrees
W = 636636;                     % lb_f
m = W/32.2;                     % slug
Jyy = 33.1e6;                   % slug-ft^2
Jxx = 18.2e6;                   % slug-ft^2
Jzz = 49.7e6;                   % slug-ft^2
Jxz = -0.35e6;                  % slug-ft^2
T_inf = 216.650;                % Kelvin
gamma = 1.4;
R = 287;                        % J/(kg-K)
a_inf = sqrt(gamma*R*T_inf);    % m/s
V_inf = M_inf*a_inf*3.28084;    % ft/s
rho = 0.000585189;              % slug/ft^3
pd = 0.5*rho*V_inf^2;           % lb_f/ft^2
S = 5600;                       % ft^2
b = 211;                        % ft
c = S/b;                        % ft
g = 32.2;                       % ft/s^2

% ---------------------
% Variables to modify
% ---------------------
hf = 45000;
tf = 50;
% ---------------------

CL = 0.521; CD = 0.0415;

global xT x0

x0 = [altitude;V_inf;0;0;0];
xT = [0;0;hf;0;0];
p0 = [1;1;1;1;1];
tspan = 0:1:tf;
options = odeset ('AbsTol',1e-6,'RelTol',1e-8);
X0 = [x0;p0];
[Tp, Xp] = ode45(@optimalTraj, tspan,X0,options);
solinit.x = Tp;
solinit.y = Xp';
solOut = bvp4c(@optimalTraj, @bcOptimalTraj, solinit);
Tp = solOut.x; Tp=Tp(:);
Xp = solOut.y; Xp = Xp';

%% Plots

subplot(3, 2, 1);
plot(Tp, Xp(:,1));
xlabel('Time (sec)');
ylabel('Altitude, Z (ft)');

subplot(3, 2, 2);
plot(Tp, Xp(:,2));
xlabel('Time (sec)');
ylabel('Horizontal Velocity, U (ft/s)');

subplot(3, 2, 3);
plot(Tp, Xp(:,3));
xlabel('Time (sec)');
ylabel('Vertical Velocity, W (ft/s)');

subplot(3, 2, 4);
plot(Tp, Xp(:,4));
xlabel('Time (sec)');
ylabel('Pitch Angle, \theta (rad)');

FT = -Xp(:,7)/(2*m);
MAC = -Xp(:,10)/(2*Jyy);

subplot(3, 2, 5);
plot(Tp, FT);
xlabel('Time (sec)');
ylabel('Thrust, F_T (lb_f)');

subplot(3, 2, 6);
plot(Tp, MAC);
xlabel('Time (sec)');
ylabel('Moment, M_{AC} (lb-ft)');

%% Functions

function [Xdot] = optimalTraj(~, X)

global m Jyy rho g CD S CL

p = X(6:10);
Z = X(1); U = X(2); W = X(3); th = X(4); Q = X(5);

alfa = atan2(W, U);
D = 0.5*rho*(U^2+W^2)*CD*S;
L = 0.5*rho*(U^2+W^2)*CL*S;

p1 = 0;
p2 = p(1);
p3 = p(2);
p4 = p(3);
p5 = p(4);
p6 = p(5);

FT = -p3/(2*m);
MAC = -p6/(2*Jyy);

Xdot(1) = -U*sin(th)+W*cos(th);
Xdot(2) = -W*Q-g*sin(th)-cos(alfa)*D/m + sin(alfa)*L/m + FT/m;
Xdot(3) = U*Q+g*cos(th)-sin(alfa)*D/m-cos(alfa)*L/m;
Xdot(4) = Q;
Xdot(5) = MAC/Jyy;

pdot = [0;...
    (2*CD*S*U^2*p3*rho - 2*Q*U*m*p4*((U^2 + W^2)/U^2)^(1/2) +...
    2*CL*S*U^2*p4*rho + CD*S*W^2*p3*rho + CL*S*W^2*p4*rho...
    - 2*U*m*p1*cos(th)*((U^2 + W^2)/U^2)^(1/2) +...
    2*U*m*p2*sin(th)*((U^2 + W^2)/U^2)^(1/2) +...
    CD*S*U*W*p4*rho - CL*S*U*W*p3*rho)/(2*U*m*((U^2 + W^2)/U^2)^(1/2));...
    (2*Q*U*m*p3*((U^2 + W^2)/U^2)^(1/2) + CD*S*U^2*p4*rho...
    - CL*S*U^2*p3*rho + 2*CD*S*W^2*p4*rho - 2*CL*S*W^2*p3*rho...
    - 2*U*m*p2*cos(th)*((U^2 + W^2)/U^2)^(1/2) -...
    2*U*m*p1*sin(th)*((U^2 + W^2)/U^2)^(1/2) + CD*S*U*W*p3*rho...
    + CL*S*U*W*p4*rho)/(2*U*m*((U^2 + W^2)/U^2)^(1/2));...
    p2*(U*cos(th) + W*sin(th)) - p1*(W*cos(th) - U*sin(th))...
    + g*p3*cos(th) + g*p4*sin(th);
    W*p3 - U*p4 - p5];

Xdot(6:10) = pdot;
Xdot = Xdot(:);

end

function bc = bcOptimalTraj(XA, XB)

global xT x0

bc(1:5) = x0(:)-XA(1:5);
bc(6:10) = [XB(4)-xT(1);XB(5)-xT(2);XB(1)-xT(3);XB(7)-xT(4);XB(3)-xT(5)];

end