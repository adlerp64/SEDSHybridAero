g = 32.17;
t = 6.5;
tStep = 0.05;
vX = 10;
vY= 700;
v = (vX^2 + vY^2)^0.5;
y = 0;
x = 0;
alpha = 0;
phi = 0;
alt = 1000;
omega = 0;
A = 3.14159 * 3 * 3;
launchFromSea = 0;

ftToM = .3048;
lbToKg = 0.4536;


dMass = 125;

xDev = 0;
%xDev = 0.02/12; %X distance of center of mass from long. axis of rocket

[T,a,P,rho] = ownAlt (alt);

%function [normalForce, cpLoc] = normalForce(alpha, rho, v)
%function [thrust, propMass, CG, m_dot] = thrustTwo(t)
%function [momInert] = inertX (propMass)
%function [dragForce] = dragForce(alpha,rho,v)
%function [dampTorq] = dampTorque (mDot, dCOM, omega)


altNew = 1000;
%results = ["t", "alt", "m", "v", "vY", "vX","aY", "aX", "alpha", "fThrust", "fDrag", "fNormal"]; 
%results = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26];
%results = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]; THIS ONE
%is big, main mark 2 one
%results = [1,2,3,4,5,6,7]
results = [1,2,3,4,5,6,7,8,9,10,11,12,13];
while altNew >= alt || alt < 1000
    alt = altNew;
    if alt < 100
        tStep = 0.01;
    elseif t > 4
            tStep = 0.01;
    else
        tStep = 0.05;
    end
    
    
    aboveSea = alt+launchFromSea;
    %[T,a,P,rho] = ownAlt (aboveSea);   INCLUDE THIS
    
    rho = 0.0408;
    %rho = rho * 32.17; %Convert slugs to lbs
    %rhoIn = rho / 12/12/12; %Convert lb/ft^3 to lb/in^3
    [fNormal, cpLoc] = normalForce ((phi-alpha),rho, v);
    %[fNormal, cpLoc] = normalForce (alpha, rho, v);
    [fThrust, propMass, dCOM, mDot] = thrustTwo (t);
    [tDamp,dampCof] = dampTorque (mDot, dCOM, omega);
    if fThrust == 0
        tDamp = 0;
    end
    
    
    %fDrag = dragForce (phi, rho, v);     INCLUDE THIS
    
    fDrag = 0.6 * 3.14159 * (3/12)^2 * rho * v^2;
    %fDrag = 0; %TEST TO SEE HOW MUCH CONTRIB

    cPDif = cpLoc - dCOM;
    cPDif = cPDif / 2;
    
    m = dMass + propMass;
    massNotForce = m;
   
    %Acceleration Equations
    fThrustX = fThrust * sin (alpha);
    fThrustY = fThrust * cos (alpha);
    
    fNormalX = fNormal * cos (alpha);
    fNormalY = fNormal * sin (alpha);
    
    tDamp = tDamp + 0;
    
    if v >= 10
        fDragX = fDrag * (vX / v);
        fDragY = fDrag * (vY / v);
    else
        fDragX = 0;
        fDragY = 0;

    end
    fX = fThrustX + fNormalX + fDragX;

    fY = fThrustY + fNormalY - fDragY - m * g;
    
    y = y + vY * tStep;
    x = x + vX * tStep;
    
    if vX == 0   
        aX = fX / massNotForce;
    else
        aX = fX / massNotForce * abs(vX) / vX;
    end
    
    aY = fY / massNotForce;   
    
    if alt < 40
        vX = 0;
        aX = 0;
    end
    
    vX = vX + aX * tStep;
    vY = vY + aY * tStep;
    v = (vX^2 + vY^2)^0.5;
    
    if t<2 && vY < 1
        vY = 0;
    end
    
    
    if vY ~= 0
        phi = atan (vX/vY);
    else
        phi = 0;
    end
    
    %Moment Equations
    dragLen = xDev * cos(alpha) + cPDif * sin (alpha);
    momInertOne = inertX(propMass);
   %function [torque, thrustTorque, dragTorque, normalForceTorque] = netTorque (alt, COPDif, xDev, alpha, phi, fDrag, fThrust, fNormal, dampTorque, omega)
    [torque,thrustTorque,dragTorque,normalForceTorque] = netTorque (alt, cPDif, xDev, alpha, phi, fDrag, fThrust, fNormal, tDamp, omega);
    
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
    %newResults = [t,alt,v,vY,vX,alphaDeg,phiDeg,aoa,omega,omegaDot,torque,thrustTorque,dragTorque,normalForceTorque, tDamp, fThrust, fDrag, fNormal,mDot,dampCof];
    %newResults = [t,alt,v,alphaDeg,cpLoc,dCOM,cpMargin];
    newResults = [t,alt,v,alphaDeg,phiDeg,aoa,omega,torque,thrustTorque,dragTorque,normalForceTorque,fDrag,fNormal];
    if t > 3
        results = [results; newResults];
    end
    
    disp (t);
    disp (v);
    disp (alt);
    disp (alphaDeg);
    %disp (alphaDeg);
    %if alt >= 40 && t > 3
        % Convert cell to a table and use first row as variable names
        %T = cell2table(c(2:end,:),'VariableNames',c(1,:))
 
        % Write the table to a CSV file
        %writetable(T,'myDataFile.csv')
        %T = table(results);
        csvwrite ('04_700fpsAt3000ftCheckingIssues.csv', results);
        %writetable (T, 'trajectoryData.txt');
    %end
    
end
      