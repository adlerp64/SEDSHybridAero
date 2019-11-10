%function [suc, transYieldLoad] = transYieldFailure (Fy, D, t, W, a, alum, ultYieldStress, margin)
%Moment of Inertia: inertX
%	Input: []
%	Output: (momInert)

%Unit conversions
% 0.0254 m = 1 in
% 0.4545 kg = 1lb

function [inertX] = momInert (propMass)

    mParafin = propMass * 3/24;
    mOx = propMass * 21/24;
    
    lbToKg = 0.4545;
    inToM = .0254;
    
    LProp = 19.65;
    propCenter = 9.5 * 12;
    LOx = 33.41;
    oxCenter = 6.5 * 12;
    
    dryCenter = 83.74;
    
    IProp = 1/12 * mParafin * lbToKg * inToM * LProp^2;                %Equations based on geometry and conversions into meters
    IProp = IProp + mParafin * lbToKg * ((propCenter - dryCenter) * inToM)^2; % Parallel axis theroem
    
    IOx = 1/12 * mOx * lbToKg * inToM * LOx^2;
    IOx = IOx + mOx * lbToKg * ((dryCenter - oxCenter) * inToM)^2;

    IDry = 207991.79 %In lb in^2
    IDry = IDry * lbToKg * inToM^2;
    
    inertX = IDry + IOx + IProp;
    
    
    
    
    
