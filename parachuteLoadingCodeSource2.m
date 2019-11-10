
%Parachute Code

W = 95;   %lbs- Dry mass 95 lbs
Cd = 1.6; %Max Cd the parachute may have, therefore largest force
Vs = 150; % ft/s - From open rocket 20 degree code
Vs = 150;
rho = 0.053172; %oz/in^3 - From apogee open rocket starting at sea level
rho = rho * 144*12/16; %oz/ft^3
rho = 0.001756; %slug/ft^3 - non open rocket
rho = rho * 32; %lb/ft^3
To = 0.09; % s - Based on elliptical graph given
So = 12.6; %ft^2 - Deployment area of parachute
g = 32.2; % ft/s^2 - Gravitational acceleration
eta = 0; %ratio of stretch area to fully deployed area
Fs = So * Cd * .5 * Vs^2 * rho/32; %units? - Force assuming drag equation at full deployment
FTerm = So * Cd * .5 * 64.96^2 * rho/32;

%Nylon Material Properties - Assumed Nylon 6 (eng. toolbox)
nylonTensMod = 2; %GPa - ranges from 2-4, took lowest
UltTensStrength = 45; % MPa- Ranges from 45-90 took lowest
yieldStrength = 45; % MPa

%Mass ratio - Uses approximation of To from graph
M = 2*W / (rho * Vs * To * Cd * So);


%Shock factor at full expansion - 1st one assumes 0 area at stretch
%Use equation 22 for if area at stretch is given
%Xo = 1/(1 + 1/(7*M))^2;
Xo = (1-eta)^2 + 2*eta * (1-eta) + eta^2;
Xo = Xo / (1 + 1/M * ((1-eta)^2/7 + (eta*(1-eta))/2 + eta^2))^2;
disp ("Xo");
disp (Xo);


%Finding Max Area Ratio
%gammaGraph2 = importdata ('r20Tr55PinBendingVariable.csv');
%graph1 = gammaGraph1.data;

maxAreaRatioGraph = importdata ('dragAreaRatioToInitElong.csv');
maxAreaRatioTable = maxAreaRatioGraph.data;
disp ("maxAreaRatioTable");
disp (maxAreaRatioTable);
disp ("Xo");
disp (Xo);
%function [k] = sortGraph (k_XAxis, k_YAxis, kGraph1, secondAxis)
areaRatio = sortGraph (Xo, 2, maxAreaRatioTable, 2);
disp ("Xo, areaRatio");
disp (Xo);
%areaRatio = 1;
disp (areaRatio);

%maxTimeRatio = Tf/To
maxTimeRatio = areaRatio^(1/6);

%velocityRatio = Vo / Vs
velocityRatio = 1 / (1 + 1/M * ((1-eta)^2/7 + eta*(1-eta)/2 + eta^2));
K = velocityRatio; %Simplification used by source

xMax = maxTimeRatio^6;
xMax = xMax / (1/K + 1/(7*M) * (maxTimeRatio^7 - 1))^2;

FMax = xMax * Fs;
xEndExtra = 16/49 * (21*M / 4)^(6/7);

disp ("M,Fs, xMax, FMax, xEndExtra, FTerm");
disp (M);
disp (Fs);
disp (xMax);
disp (FMax);
disp (xEndExtra);
disp (FTerm);


%Eo = elongation at To

%Eo = Xo * Fs / Fc * Emax






