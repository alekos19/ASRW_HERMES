%% PARAMETERS

vesselCapacity  = 2.0;
waterCapacity   = 1.0;
vesselDiameter  = 95;
nozzleDiameter  = 10;
rocketDryMass   = 0.20;
pressureGauge   = 6.0;

param = struct();
param.k                 = 1.40;
param.g                 = 9.81;
param.pa                = 1e5;
param.ro                = 1000;
param.rho               = 1.225;
param.R                 = 286.9;
param.pI                = pressureGauge * 1e5;
param.tempI             = 293;
param.waterLossFactor   = 0.20;
param.airLossFactor     = 0.30;
param.C                 = param.pI * ((vesselCapacity - waterCapacity) / 1000)^param.k;
param.Ar                = pi * (vesselDiameter  / 1000)^2 / 4;    
param.An                = pi * (nozzleDiameter  / 1000)^2 / 4;                       
param.m0                = rocketDryMass;
param.fi                = deg2rad(90.0);

% define thrust state
Vw     = waterCapacity / 1000;
Va     = (vesselCapacity - waterCapacity) / 1000; 
p      = param.pI;
temp   = param.tempI;
mA     = (p + param.pa) / (param.R * temp) * Va; 
tank = [Vw Va p temp mA];

function [T, tank_dot] = thrust(param, tank)
    Vw     = tank(1);
    Va     = tank(2); 
    p      = tank(3);
    temp   = tank(4);
    mA     = tank(5);

    if Vw > 0
        ve = sqrt(2*p*(1-param.waterLossFactor) / param.ro);
        Vw_dot = -param.An * ve; 
        Va_dot = -Vw_dot;

        p_dot = -param.k * param.C * Va^(-param.k-1) * Va_dot;
        % optional hydrostatic p_dot = p_dot + param.ro * (param.g + a)/param.Ar * Vw_dot

        n = (param.k-1)/param.k;
        temp_dot = param.tempI * n * ((p + param.pa)/(param.pI + param.pa))^(n-1) * (1/(param.pI + param.pa))*p_dot;

        m_dot = Vw_dot * param.ro;
        T = ve*(-m_dot)+p*(1-param.waterLossFactor)*param.An;
        tank_dot = [Vw_dot Va_dot p_dot temp_dot 0];
    elseif p > 0
        roApE  = mA/Va;
        ve     = sqrt(temp * param.R * (2*param.k/(param.k-1)) * (1 - (param.pa/(p+param.pa))^((param.k-1)/param.k)));
        mA_dot = -ve * param.An * roApE;

        p_dot = -param.k * param.C * Va^(-param.k-1) * (-mA_dot / roApE);
        n = (param.k - 1)/param.k;
        temp_dot = param.tempI * n * ((p + param.pa)/(param.pI + param.pa))^(n-1) * (1/(param.pI + param.pa)) * p_dot;

        T  = ve * (1-param.airLossFactor) * (-mA_dot) + p * param.An;
        tank_dot = [0 0 p_dot temp_dot mA_dot];
    else
        T = 0;
        tank_dot = [0 0 0 0 0];
    end
end

function [X_dot,tank_dot, T] = dynamics(X, tank, fi, param)
    x = X(1);
    x_dot = X(2);
    y = X(3);
    y_dot = X(4);
    
    m = mass(tank, param);
    w = m * param.g;
    
    [T, tank_dot] = thrust(param, tank);
    y_dotdot = (T*sin(fi) - w) / m;
    x_dotdot = (T*cos(fi)) / m;
    
    X_dot = [x_dot x_dotdot y_dot y_dotdot];
end

function m = mass(tank, param)
    m = param.m0 + tank(1) * param.ro + tank(5);
end

%% SIMULATION

dt      = 0.001;
t_max   = 7;
t_vec   = 0:dt:t_max;
N       = length(t_vec);

thrust_log      = zeros(N,1);
pressure_log    = zeros(N,1);
X_log = zeros(N,4);
a_log = zeros(N,2);
i_cuttoff = -1;

X = [ 0 0 0 0 ];
fi = deg2rad(90.0);
for i=1:(N-1)
    [kx1,kt1]   = dynamics(X, tank, fi, param);
    [kx2,kt2]   = dynamics(X+0.5*dt*kx1, tank+0.5*dt*kt1, fi, param);
    [kx3,kt3]   = dynamics(X+0.5*dt*kx2, tank+0.5*dt*kt2, fi, param);
    [kx4,kt4,T] = dynamics(X+1*dt*kx3, tank+1*dt*kt3, fi, param);

    X = X + (dt/6) * (kx1 + 2*kx2 + 2*kx3 + kx4);
    tank = tank + (dt/6) * (kt1 + 2*kt2 + 2*kt3 + kt4);

    if T == 0 && i_cuttoff < 0
        i_cuttoff = i;
    end

    thrust_log(i)   = T;
    pressure_log(i) = tank(3) / 1e5;

    X_log(i,:) = X;
    a_log(i,:) = [kx4(2) / param.g kx4(4) / param.g];
end

subplot(3,2,1);
plot(t_vec(1:i_cuttoff), thrust_log(1:i_cuttoff));
title("Thrust (N)");
grid on

subplot(3,2,2);
plot(t_vec(1:i_cuttoff), pressure_log(1:i_cuttoff));
title("Pressure (bar)");
grid on

subplot(3,2,3);
plot(t_vec, X_log(:,4));
title("Vy (m/s)");
grid on

subplot(3,2,4);
plot(t_vec, X_log(:,3));
title("Y (m)");
grid on

subplot(3,2,5);
plot(t_vec(1:i_cuttoff), a_log((1:i_cuttoff),2));
title("ay (g)");
grid on