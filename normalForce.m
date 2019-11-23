function [normalForceExp, cpLoc] = normalForce(alpha, rho, v)


%% Init
% Nose Cone Parameters
l_n = 2;                                % [ft] length of nose cone
d_n = 0.5;                              % [ft] diameter of nose cone

% Body Parameters
l_b = 140/12;                               % [ft] length of body
d_b = 0.5;                              % [ft] diameter of body
x_b = 2;                                % [ft] location of body from tip
l_tr = l_n + l_b;                       % [ft] length of whole rocket

% Fin Parameters
n = 3;                                  % [ ]  number of fins
l_r = 10/12;                            % [ft] fin root chord length
l_t = .77/12;                           % [ft] fin tip chord length
l_w = 6.781/12;                         % [ft] fin sweep length
l_s = 7.75/12;                          % [ft] fin span length
x_f = 154/12;                           % [ft] location of fins from tip
d_f = 0.5;                              % [ft] body diameter at fins
l_m = sqrt(l_s^2 + (l_r-l_t)^2);        % [ft] fin mid-chord length
l_ts = 2*l_s + d_f;                     % [ft] total span of fins

A_ref = 1/4*pi()*d_n^2;                 % [ft] reference area
A_p = d_b*l_b;                          % [ft*ft] planform area of body

%% Normal Foce Calc:

% Stability derivative for nose cone
% Assumption: either conical, ogive, or parabolic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEED TO FIND SOURCE FOR HAACK SERIES%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_n_a_n = 2;

% Coefficient for fin-body interference
K_fb = 1 + (d_f/2)/(l_s+d_f/2);
% Stability derivative for fins
c_n_a_f = K_fb*(4*n*(l_s/d_n)^2)/(1+sqrt(1+(2*l_m/(l_r+l_t))^2));

% Pseudo stability derivative for body-lift
K = 1;  % Needs to be between 1 and 1.5, Galejs [1999] uses 1 with good results
c_n_a2 = @(alpha) K*A_p/A_ref*alpha;
%disp ("Stability derivative for bodyLift");
normalForceBodyLift = 0.5 * rho * v^2 * A_ref * c_n_a2(alpha) * alpha;
%disp (normalForceBodyLift);

% Total stability derivative
c_n_a = @(alpha) c_n_a_n + c_n_a_f + c_n_a2(alpha);
c_n_a = @(alpha) c_n_a_n + c_n_a_f; %No body force

% Normal force coefficient
c_n = @(alpha) c_n_a(alpha)*alpha;

% Normal Force
normalForceExp = 0.5 * rho * v^2 * A_ref * c_n(alpha);

%% CP location Calc:

% Center of pressure location for nose cone
x_cp_n = 0.5*l_n;                       % [ft] cp location for von karman nose cone

% Center of pressure location for fins
x_cp_f = x_f + (l_m*(l_r+2*l_t))/(3*(l_r+l_t)) + 1/6*(l_r+l_t-(l_r*l_t)/(l_r+l_t));

% Center of pressure location for body lift (center of planform area)
x_cp_b = 0.5*(0.40355 + (x_b + 0.5*l_b));

% Center of pressure of rocket
x_cp = @(alpha) (c_n_a_n*x_cp_n + c_n_a_f*x_cp_f + c_n_a2(alpha)*x_cp_b)/c_n_a(alpha);

cpLoc = x_cp(alpha);                    % [ft] location of center of pressure from tip

%disp("fins");
normalForceFins = 0.5 * rho * v^2 * A_ref * c_n_a_f * alpha;
%disp (normalForceFins);
%disp ("nose");
normalForceNose = 0.5 * rho * v^2 * A_ref * 2 * alpha;
%disp (normalForceNose);
%disp ("body")
%disp (normalForceBodyLift);
%disp ("total")
%disp (normalForceExp);







end