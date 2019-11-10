function [thrust, propMass, CG, m_dot] = thrust(t)

%% Init
m_prop = 27;                % [lbm] propellant mass % Updated value is 24.
m_dry = 125;                % [lbm] dry mass of rocket 
x_dry = 5.082569431;        % [lbm] dry CG location
x_tank = 6.5;               % [ft] location of middle of oxidizer tank
x_fg = 9.5;                 % [ft] location of middle of fuel grain

Isp = 165;                  % [s] Specific Impulse
mdot_inital = 9.39;         % [lbm/s] Inital Mass Flow Rate
g = 32.2;                   % [ft/s/s] Newton Gravitational Constant
mdot_average = 4.21;        % [lbm/s] Average Mass Flow Rate
%% Calculations

% Thrust at time t (interpolation from thrust curve)
thrust_table = importdata('hybrid_thrust_curve.eng',' ');
T = @(t) interp1(thrust_table(:,1), thrust_table(:,2),t,'linear', 0);
thrust = T(t) / 4.4482;     % [lbf] conversion from N to lbf

% Propellant mass
% Assumption: change in mass is proportional to change in impulse
propMass = @(t) m_prop - m_prop * integral(T, 0, t) / integral(T, 0, Inf);


% Mass flow rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need better numbers from propulsion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% m_dot = @(t) T(t) / (Isp * g);
% m_dot = @(t) interp1([0,thrust_table(end,1)], [mdot_inital, 0], t, 'linear',0);
m_dot = mdot_average;

% Center of Gravity
% Assumption: propellants are point masses and rate of change is equal
x_prop = 0.5*(x_tank + x_fg);

% CG equation
CG = (m_dry*x_dry + propMass(t)*x_prop) / (m_dry + m_prop);

%% Outputs

propMass = propMass(t);

end