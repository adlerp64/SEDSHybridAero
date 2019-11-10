
g = 32.17;
t = 0;
tStep = 0.01;
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
A = 3.14159 * 3 * 3;
launchFromSea = 0;
launchRodAngle = 5 * 0.01745;

ftToM = .3048;
lbToKg = 0.4536;


dMass = 125;

%xDev = 0;
xDev = 0.02/12; %X distance of center of mass from long. axis of rocket

[T,a,P,rho] = ownAlt (alt);

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
results = [1,2,3,4,5,6,7,8,9];
while altNew >= alt || alt < 1000
    alt = altNew;
    if alt < 100
        tStep = 0.01;
    elseif fThrust ==0
            tStep = 0.001;
            %tStep = 0.01;
    else
        tStep = 0.05;
    end
    if alt <55
        if alt > 45
            %vX = 5;
            vX = vX * 1;
            
        end
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
    %function [dDampT] = dragDamp (omega, rho, com, l)
    dDampT = dragDamp(omega,rho, dCOM, 10);
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

    %if alphaDeg > 2
    %   g =0;
    %end
    fY = fThrustY + fNormalY - fDragY - m * g;
    
    y = y + vY * tStep;
    x = x + vX * tStep;
    
    if vX == 0   
        aX = fX / massNotForce;
    else
        aX = fX / massNotForce * abs(vX) / vX;
    end
    
    aY = fY / massNotForce;   
    
    if alt < 30*cos(launchRodAngle)
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
    %torque  = torque + dDampT;
    
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
    newResults = [t,alt,v,vY,vX,alphaDeg,aoa,x,y];
    if t > 3
        results = [results; newResults];
    end
    
    disp (t);
    disp (v);
    disp (alt);
    disp (alphaDeg);
    disp (aoa);
    %disp (alphaDeg);
    if alt >= 40 && t > 3
        % Convert cell to a table and use first row as variable names
        %T = cell2table(c(2:end,:),'VariableNames',c(1,:))
 
        % Write the table to a CSV file
        %writetable(T,'myDataFile  .csv')
        %T = table(results);
        csvwrite ('05_LaunchRodAdded.csv', results);
        %writetable (T, 'trajectoryData.txt');
    end
    
end
      