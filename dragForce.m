function [drag] = dragForce(alpha,rho,v)
%DRAGFORCE computation of drag force using Barrowman Equations
%   dragForce:      [lbf] force of drag
%   alpha:          [rad] angle of attack
%   rho:            [lbm/ft/ft/ft] air density
%   v:              [ft/sec] relevent velocity

%% Init
% Nose Cone Parameters
l_n = 2;                                % [ft] length of nose cone
d_n = 0.5;                              % [ft] diameter of nose cone

% Body Parameters
l_b = 10;                               % [ft] length of body
d_b = 0.5;                              % [ft] diameter of body
x_b = 2;                                % [ft] location of body from tip
l_tr = l_n + l_b;                       % [ft] length of whole rocket

% Fin Parameters
n = 3;                                  % [ ]  number of fins
l_r = 10/12;                            % [ft] fin root chord length
l_t = 3.172/12;                         % [ft] fin tip chord length
l_w = 6.781/12;                         % [ft] fin sweep length
l_s = 5.512/12;                         % [ft] fin span length
x_f = 154/12;                           % [ft] location of fins from tip
d_f = 0.5;                              % [ft] body diameter at fins
l_m = sqrt(l_s^2 + (l_r-l_t)^2);        % [ft] fin mid-chord length
l_ts = 2*l_s + d_f;                     % [ft] total span of fins

A_ref = 1/4*pi()*d_n^2;                 % [ft] reference area

%% Calculations

% Zero angle of attack drag coefficient
c_d_0 = 0.6;    % Assumption (pretty accurate)

delta_table = csvread('delta.csv',1);   % Obtained from Mandell et al, 1973
eta_table = csvread('eta.csv',1);       % Obtained from Mandell et al, 1973

% Interpolation of table data
delta = @(alpha) interp1(delta_table(:,1),deg2rad(delta_table(:,2)),alpha,'linear',0);
eta = @(alpha) interp1(eta_table(:,1),deg2rad(eta_table(:,2)),alpha,'linear',0);

% Additional drag at angle of attack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS IS BEING REALLY CRAPPY also domain is only from 4 to 20 deg%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_d_b_a = @(alpha) 2*delta(alpha)*alpha^2 + (3.6*eta(alpha)*(1.36*l_tr - 0.55*l_n))/(pi()*d_b)*alpha^3;

% Additional data on fins
R_s = l_ts/d_f;                                 % [ ] fin section ratio
k_fp = 0.8065*R_s^2 + 1.1553*R_s;               % [ ] fin-body interference coefficient
k_bf = 0.1935*R_s^2 + 0.8174*R_s + 1;           % [ ] body-fin interference coefficient
A_fe = n*0.5*(l_r + l_t)*l_s;                   % [ft*ft] total exposed area of fin
A_fp = n*A_fe + 0.5*d_f*l_r;                    % [ft*ft] total planform area of fin

% Coefficient of alpha drag on fins
c_d_f_a = @(alpha) alpha^2*(1.2*A_fp*4/(pi()*d_f^2) + 3.12*(k_fp+k_bf-1)*(A_fe*4/(pi()*d_f^2)));

% Coefficient of drag as a function of alpha
c_d = @(alpha) c_d_0 + c_d_b_a(alpha) + c_d_f_a(alpha);

% Drag force as a function of alpha
%drag = rho*v^2*A_ref*c_d(alpha);       % [lbf] force due to drag
drag = rho*v^2*A_ref*c_d_0;
drag = abs(drag);
%dragForce = c_d_0 * rho*v^2*A_ref*0.5*.25; %Replacement dragForce for Cd = 0.5 multiplied by 0.5 again for error
%disp("Drag Force");
%disp (dragForce);

end