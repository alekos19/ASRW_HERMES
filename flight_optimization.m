STEPS=500;

cf_rise = 0.08;
cf_fall = 0.8;
Sref = .1 ;    
vesselCapacity  = 2.0;
nozzleDiameter  = 9;
waterLossFactor = 0.20;
airLossFactor   = 0.30;
k = 1.40;
g = 9.81;
pa = 100000;
ro = 1000;
R  = 286.9;
temp = 293;
    
height_log = zeros(STEPS,1);
time_log   = zeros(STEPS,1);

for pressureGauge=5:10
    for RUN=1:STEPS
        fprintf("Running pressure="+pressureGauge+" "+RUN+"/"+STEPS+"\n");
    
        waterCapacity = vesselCapacity * (RUN/STEPS);
        rocketDryMass   = 2;
        tmax = 20;
        dt   = 0.001;
        N = tmax / dt;
        [T_log, ve_log, p_log, m_log] = thrust(vesselCapacity, waterCapacity, nozzleDiameter, rocketDryMass, pressureGauge, waterLossFactor, airLossFactor, k, pa, ro, R, temp, N, dt);
        
        % Κατάσταση: [y; vy] — μόνο κατακόρυφος άξονας
        X_log = zeros(2, N);
        X_log(:,1) = [0; 0];  % y=0, vy=0
        
        for i = 1:(N-1)
            m   = m_log(i);
            T   = T_log(i);
            rho = 1.225;
            Xi  = X_log(:,i);
        
            k1 = dt * dynamics(Xi,         T, m, g, cf_rise, cf_fall, rho, Sref);
            k2 = dt * dynamics(Xi+0.5*k1,  T, m, g, cf_rise, cf_fall, rho, Sref);
            k3 = dt * dynamics(Xi+0.5*k2,  T, m, g, cf_rise, cf_fall, rho, Sref);
            k4 = dt * dynamics(Xi+k3,      T, m, g, cf_rise, cf_fall, rho, Sref);
        
            X_log(:,i+1) = Xi + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
        
            if X_log(1, i+1) < 0
                height_log(RUN) = max(X_log(1,:));
                time_log(RUN) = i * dt;
                break
            end
        end
    end

    figure(pressureGauge);
    sgtitle(pressureGauge);
    x = 1:STEPS;

    [maxtime,idx]=max(time_log);
    subplot(3,1,1);
    plot(x, time_log);
    title("TIME MAX="+max(time_log)+", AT="+(idx/STEPS));
    xlabel("RUN");
    ylabel("DURATION (s)");

    [maxheight,idx]=max(height_log);
    subplot(3,1,2);
    plot(x, height_log);
    title("HEIGHT MAX="+maxheight+", AT="+(idx/STEPS));
    xlabel("RUN");
    ylabel("HEIGHT (m)");

    score = time_log + height_log;
    [maxscore,idx] = max(score);
    subplot(3,1,3);
    plot(x,score);
    title("SCORE MAX="+maxscore+", AT="+(idx/STEPS));
    xlabel("RUN");
    ylabel("SOCRE");
end

% --- Dynamics: μόνο κατακόρυφος άξονας ---
function dXdt = dynamics(X, T, m, g, cf_rise, cf_fall, rho, Sref)
    vy       = X(2);
    if vy >= 0
        cf = cf_rise;
    else
        cf = cf_fall;
    end

    w        = m * g;
    F        = cf * 0.5 * rho * Sref * vy^2 * sign(vy); 
    y_dotdot = (T - F  - w) / m;
    dXdt     = [vy; y_dotdot];
end



function [T_log, ve_log, p_log, m_log] = thrust(vesselCapacity, waterCapacity, nozzleDiameter, rocketDryMass, pressureGauge, waterLossFactor, airLossFactor, k, pa, ro, R, temp, N,dt)
    V   = (vesselCapacity - waterCapacity) / 1000; % Αρχικός όγκος αέρα (m³)
    Vw  = waterCapacity / 1000;                    % Αρχικός όγκος νερού (m³)
    An  = pi * (nozzleDiameter  / 1000)^2 / 4;     % Διατομή ακροφυσίου (m²)
    p   = pressureGauge * 1e5;                     % Gauge πίεση (Pa)
    pI    = p;                                     % Αρχική πίεση αναφοράς (Pa)
    tempI = temp;                                  % Αρχική θερμοκρασία αναφοράς (K)
    mA    = (p + pa) / (R * temp) * V;             % Μάζα συμπιεσμένου αέρα (kg)
    m0    = rocketDryMass;                         % Ξηρή μάζα (kg)
    m     = m0 + Vw * ro + mA;                     % Συνολική αρχική μάζα (kg)
    C     = p * V^k;                               % Αδιαβατική σταθερά C = p·V^γ
    
    T_log = zeros(N, 1);
    ve_log= zeros(N, 1);
    p_log = zeros(N, 1);
    m_log = ones(N,1)*m0;
    
    idx = 0;  
    while idx < N && m > m0 + mA
        ve = sqrt(2 * p * (1 - waterLossFactor) / ro); % Exhaust velocity νερού (m/s)
        dm = An * ro * ve * dt;                         % Μάζα νερού που εξέρχεται (kg)
    
        m  = m  - dm;          % Μείωση συνολικής μάζας
        V  = V  + dm / ro;     % Αύξηση V_air
        Vw = Vw - dm / ro;     % Μείωση V_water
        
        p    = C * V^(-k);
        temp = tempI * ((p + pa) / (pI + pa))^((k-1)/k);
    
        T = ve * (dm/dt) + p * (1 - waterLossFactor) * An; % Ώση (N)
        
        idx = idx + 1;
        m_log(idx)    = m;
        T_log(idx)    = T;
        ve_log(idx)   = ve;
        p_log(idx)    = (p + pa) / 1e5;
    end
    
    roApE = mA / V;      % Αρχική εκτίμηση πυκνότητας αέρα στο λαιμό
    while idx < N && mA > 0
        % --- Exhaust velocity: ισεντροπική εκτόνωση από p_chamber → pa ---
        ve = sqrt(temp * R * (2*k/(k-1)) * (1-(pa/(p+pa))^((k-1)/k)));
    
        % --- Θερμοκρασία & ταχύτητα ήχου στο λαιμό (critical conditions) ---
        tempE  = temp * 2 / (k + 1);       % T_throat = T · 2/(γ+1) (isentropic)
        vSound = sqrt(k * R * tempE);      % a_sound = √(γ·R·T_throat)
        M      = ve / vSound;              % Mach στο λαιμό
    
        if M >= 1
            % Choked flow: ροή φτάνει ταχύτητα ήχου στο λαιμό (NASA formula)
            dm = An * p / sqrt(temp) * sqrt(k/R) * ((k+1)/2)^(-(k+1)/(2*(k-1))) * dt;
            ve = vSound;              % ve δεν μπορεί να υπερβεί vSound ισεντροπικά
        else
            % Unchoked flow: κλασική συνέχεια
            dm = ve * An * roApE * dt;
        end
    
        roApE = dm / (ve * dt * An);  % Ενημέρωση πυκνότητας αέρα στο λαιμό (kg/m³)
    
        % --- Ενημέρωση μάζας ---
        m  = m  - dm;                 % Μείωση συνολικής μάζας
        mA = mA - dm;                 % Μείωση μάζας αέρα
    
        % --- Πίεση: εικονική μεταβολή όγκου για να εξαχθεί νέα πίεση ---
        C  = p * V^k;                 % Ενημέρωση αδιαβατικής σταθεράς
        V  = V + dm / roApE;          % Εικονική αύξηση V (αέρας που εξέρχεται)
        p  = C * V^(-k);              % Νέα πίεση μέσω pV^γ = C
        temp = tempI * ((p + pa) / (pI + pa))^((k-1)/k); % Αδιαβατική θερμοκρασία
        V  = V - dm / roApE;          % Επαναφορά: το δοχείο έχει ΣΤΑΘΕΡΟ όγκο!
    
        % --- Effective ve με loss factor (μέρος ροής δεν βγαίνει κατακόρυφα) ---
        ve_eff = ve * (1 - airLossFactor);
    
        % --- Ώση ---
        T = ve_eff * (dm/dt) + p * An;  % [momentum] + [pressure thrust]
        
        idx = idx + 1;
        m_log(idx)    = m;
        T_log(idx)    = T;
        ve_log(idx)   = ve;
        p_log(idx)    = (p + pa) / 1e5;
    end
end


