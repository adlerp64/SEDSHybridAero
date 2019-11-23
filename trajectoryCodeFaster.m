clc; clear all;
g = 32.17;
resultsCount = 0;
t = 0;
tStep = 0;
tStepMost = 0.0001;
vX = 0;
vY= 0;
v = 0;
y = 0;
x = 0;
alpha = 0;
alphaDeg = alpha / 0.01745; 
phi = 0;
alt = 0;
omega = 0;
wind = 5; %ft/s assuming constant for now
A = 3.14159 * 3 * 3/144;
launchFromSea = 0; %Altitude relative to sea level of starting location
launchRodAngle = 5 * 0.01745; %Launch rod angle deg * conversion to rad
launchRodLength = 30;
%launchRodAngle = 0;
ftToM = .3048;
lbToKg = 0.4536;


dMass = 95;

xDev = 0;
%xDev = 0.02/12; %X distance of center of mass from long. axis of rocket

[T,a,P,rho] = ownAlt (alt);

%thrustTable
%altData = importdata ('altitudeTable.csv');
%altTable = altData.data;
thrustTable1 = importdata ('ITARThrust.csv');
thrustTable = thrustTable1.data;

%function [normalForce, cpLoc] = normalForce(alpha, rho, v)
%function [thrust, propMass, CG, m_dot] = thrustTwo(t)
%function [momInert] = inertX (propMass)
%function [dragForce] = dragForce(alpha,rho,v)
%function [dampTorq] = dampTorque (mDot, dCOM, omega)


altNew = 0;
%results = ["t", "alt", "m", "v", "vY", "vX","aY", "aX", "alpha", "fThrust", "fDrag", "fNormal"]; 
%results = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26];
%results = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];
%results = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
%results = [1,2,3,4,5,6,7,8,9];
results = zeros (100000, 20);
while altNew >= alt || alt < 1000
%while altNew >= alt && t < 11
    alt = altNew;
    tStep = tStepMost;
    disp ("t in main body");
    disp (t);

    aboveSea = alt+launchFromSea;
    %[T,a,P,rho] = ownAlt (aboveSea);   INCLUDE THIS
    
    rho = 0.0408;
    %rho = rho * 32.17; %Convert slugs to lbs
    %rhoIn = rho / 12/12/12; %Convert lb/ft^3 to lb/in^3
    [fNormal, cpLoc] = normalForce ((phi-alpha),rho, v);
    %[fNormal, cpLoc] = normalForce (alpha, rho, v);
    %[fThrust, propMass, dCOM, mDot] = thrustTwo (t);
    %function [thrust, propMass, CG, m_dot] = thrustValues (t, thrustTable)
    [fThrust,propMass,dCOM,mDot] = thrustValues (t,thrustTable);
    [tDamp,dampCof] = dampTorque (mDot, dCOM, omega);
    %function [dDampT] = dragDamp (omega, rho, com, l)
    dDampT = dragDamp(omega,rho, dCOM, 10);
    if fThrust == 0
        tDamp = 0;
    end
    
    
    %fDrag = dragForce (phi, rho, v);     INCLUDE THIS
    %drag = Cd * 1/2 rho v^2 A
    fDrag = 0.6 *0.5 * 3.14159 * (3/12)^2 * rho * v^2;
    %fDrag = 0; %TEST TO SEE HOW MUCH CONTRIB

    cPDif = cpLoc - dCOM;
    %cPDif = cPDif / 2; %Why is this here? Why divide the cPDif by 2?
    
    m = dMass + propMass;
    massNotForce = m;
    
    %Multiplying all forces by g
    
    %fThrust = fThrust * g;
    %fDrag = fDrag * g;
    
    
    
   
    %Acceleration Equations
    fThrustX = fThrust * sin (alpha);
    fThrustY = fThrust * cos (alpha);
    
    
    %REPALCE THIS
    %fNormalX = fNormal * cos (alpha);
    %fNormalY = fNormal * sin (alpha);
    fNormalX = 0;
    fNormalY = 0;
    
    tDamp = tDamp + 0;
    
    if v >= 10
        fDragX = fDrag * (vX / v);
        fDragY = fDrag * (vY / v);
    else
        fDragX = 0;
        fDragY = 0;

    end
    fX = (fThrustX*g + fNormalX*g - fDragX);

    %if alphaDeg > 2
    %   g =0;
    %end
    
    % a = -g - Cd*.5*getRho(h)*v^2*(pi*rocket_radius^2)/(m2)
    %g*Isp*log(m2/m1)
    % v = v + g*Isp*log(m2/m1) - g*dt - Cd*.5*getRho(h)*v^2*(pi*rocket_radius^2)/(m2)*dt;
    
    fY = (fThrustY*g + fNormalY*g - fDragY - m*g);
    
    y = y + vY * tStep;
    x = x + vX * tStep;
    
    if vX == 0   
        aX = fX / massNotForce;
    else
        aX = fX / massNotForce * abs(vX) / vX;
    end
    
    aY = fY / massNotForce;
    
    if alt < launchRodLength*cos(launchRodAngle)
        alpha = launchRodAngle;
        omega = 0;
        aTotal = (aY^2 + aX^2)^0.5;
        aYLaunchMax = aX / tan(launchRodAngle);
        aXLaunchMax = tan(launchRodAngle) * aY;
        if aY > aYLaunchMax
            aY = aYLaunchMax;
        elseif aX > aXLaunchMax
                aX = aXLaunchMax;
        end
    end
    
    vX = vX + aX * tStep;
    vY = vY + aY * tStep;
    v = (vX^2 + vY^2)^0.5;
   
    if vY ~= 0
        phi = atan ((vX-wind)/vY);
    else
        phi = 0;
    end
    
    %Moment Equations
    dragLen = xDev * cos(alpha) + cPDif * sin (alpha);
    momInertOne = inertX(propMass);
    %fNormal = 0;
   %function [torque, thrustTorque, dragTorque, normalForceTorque] = netTorque (alt, COPDif, xDev, alpha, phi, fDrag, fThrust, fNormal, dampTorque, omega)
    
   [torque,thrustTorque,dragTorque,normalForceTorque] = netTorque (alt, cPDif, xDev, alpha, phi, fDrag, fThrust*g, fNormal, tDamp, omega);
   % [torque,thrustTorque,dragTorque,normalForceTorque] = netTorque (alt, cPDif, xDev, alpha, phi, 0, 0, fNormal*g, 0, omega);
   %torque  = torque + dDampT;
    %torque = 0;
    
    alpha = alpha + omega * tStep; %Does this so that the alpha is based on the current omega, next alpha is impacted by acceleration
    
    omegaDot = torque / momInertOne;
    omega = omega + omegaDot * tStep;
    
    altNew = alt + vY * tStep;
    
    t = t + tStep;
    
    alphaDeg = alpha / 0.01745; 
    phiDeg = phi / 0.01745;

    
    cpMargin = (cpLoc - dCOM)/0.5;
    dynamicPressure = .5 * rho * v^2;
    aoa = alphaDeg - phiDeg;
    %newResults = [t, alt, m, v, vY, vX, aY, aX, alphaDeg, fThrust, fDrag, fNormal, fY, aY, fX, aX, torque, omegaDot, omega, tDamp, cPDif, momInertOne, propMass, rho, dynamicPressure,]; 
    %newResults = [t,alt,v,vY,vX,alphaDeg,phiDeg,aoa,omega,omegaDot,torque,thrustTorque,dragTorque,normalForceTorque, tDamp, dDampT, fThrust, fDrag, fNormal,mDot,dampCof,momInertOne];
    %newResults = [t,alt,v,vY,vX,alphaDeg,aoa,fDrag, fNormal,dragTorque,normalForceTorque, torque, omega, cPDif,g];
    %newResults = [t,alt,v,vY,vX,alphaDeg,aoa,x,y,cPDif];
    %propMass,dCOM,mDot
    newResults = [t,alt,v,vY,alphaDeg,aoa,fThrust,torque,thrustTorque,dragTorque,tDamp,normalForceTorque,omega,fY,aY,phi,vX,fX,aX,x];
    if t > 0.1
        results (resultsCount,:) = newResults;
    end
%      disp ("t");
%      disp (t);
%      disp (v);
%      disp (alt);
%      disp (alphaDeg);
%      disp (aoa);

        
    if t > 5-tStep && t < 5 + tStep
        csvwrite ('19_smallTStepAgain', results);
    end
    if t > 10-tStep && t < 10+tStep
        csvwrite ('19_smallTStepAgain', results);
    end
     if t > 15-tStep && t < 15+tStep
        csvwrite ('19_smallTStepAgain', results);
     end
    if t > 20-tStep && t < 20.00+tStep
        csvwrite ('19_smallTStepAgain.csv', results);
     end
    resultsCount = resultsCount + 1;
    
end
csvwrite ('19_smallTStepAgain.csv', results);

endCell = t * 1/tStep;
instantBefore = round(endCell);
instantBefore = instantBefore - 10;
pointOneBefore = instantBefore - 90;



figure(1)

subplot (2,1,1)
plot (results(1:instantBefore,1), results(1:instantBefore,2));
xlabel ('Time(s)');
ylabel ('alt (ft)');
title ('Altitude Graph');

subplot (2,1,2)
plot (results(1:instantBefore,1), results(1:instantBefore,3));
xlabel ('Time(s)');
ylabel ('v');
title ('v Graph');
 
%newResults = [t,alt,v,vY,alphaDeg,aoa,fThrust,torque,thrustTorque,dragTorque,tDamp,normalForceTorque,omega,fY,aY,phi];

figure (2)
subplot (2,1,1)
plot (results(1:pointOneBefore,1), results(1:pointOneBefore,5));
xlabel ('Time(s)');
ylabel ('alpha deg');
title ('alpha');

subplot (2,1,2)
plot (results(1:pointOneBefore,1), results(1:pointOneBefore,6));
xlabel ('Time(s)');
ylabel ('aoa');
title ('aoa');

figure (3)
subplot (2,1,1)
plot (results(1:instantBefore,1), results(1:instantBefore,18));
xlabel ('Time(s)');
ylabel ('fX');
title ('fX');


subplot (2,1,2)
plot (results(1:instantBefore,1), results(1:instantBefore,14));
xlabel ('Time(s)');
ylabel ('fY');
title ('fY');

figure (4)
subplot (2,1,1)
plot (results(1:instantBefore,1), results(1:instantBefore,4));
xlabel ('Time(s)');
ylabel ('vY');
title ('vY');

subplot (2,1,2)
plot (results(1:instantBefore,1), results(1:instantBefore,17));
xlabel ('Time(s)');
ylabel ('vX');
title ('vX');

figure (5)
subplot (2,1,1)
plot (results(1:instantBefore,20), results(1:instantBefore,2));
xlabel ('X position');
ylabel ('altitude');
title ('Position Graph');

plot (results(1:instantBefore,1), results(1:instantBefore,18));
xlabel ('Time(s)');
ylabel ('fX');
title ('fX');


% subplot (2,1,2)
% plot (results(1:instantBefore,1), results(1:instantBefore,22));
% xlabel ('Time(s)');
% ylabel ('drag');
% title ('drag force graph');
% subplot(5,1,1)
% plot(t_t, h_t);
% xlabel('Time (s)');
% ylabel('Height (ft)');
% title('Height(t)');
% 
% subplot(5,1,2)
% plot(t_t, v_t);
% xlabel('Time (s)');
% ylabel('Velocity (ft/s)');
% title('Velocity(t)');
% 




      