
function [torque, thrustTorque, dragTorque, normalForceTorque] = netTorque (alt, COPDif, xDev, alpha, phi, fDrag, fThrust, fNormal, dampTorque, omega)

momLength = (COPDif^2 + xDev^2)^0.5;
momAngle = alpha - atan(xDev/COPDif); %-aTan because counterclockwise is positive

%function [torque] = torqueFromForce (angleForce, magnitudeForce, angleLength, magnitudeLength)

if fThrust == 0
    dampTorque = 0;
end

dragTorque = -1*torqueFromForce ((phi),fDrag, momAngle, momLength);
%dragTorque = 0;
thrustTorque = torqueFromForce (alpha, fThrust, momAngle, momLength);
normalForceTorque = torqueFromForce ((3.14159/2), fNormal, 0, COPDif);
if (phi-alpha) == 0
    normalForceTorque = normalForceTorque * 1;
else
    %normalForceTorque = normalForceTorque *-1* abs((phi-alpha)) / (phi-alpha);
end

if alt < 40  %launch rail length
        torque = 0;
        
else
    if omega == 0 
        dampTorque = 0;
    end
    %normalForceTorque = 0;
    torque = thrustTorque + dragTorque + normalForceTorque + dampTorque;
            %torque = fNormal * cPDif * -1*(abs(alpha)/alpha) + tDamp * (abs(omega) / omega) + fThrust * xDev * (abs(alpha) / (alpha))+ fDrag * dragLen * (abs(alpha) / (alpha));
end
disp ("momAngle");
disp (momAngle/0.01745);

end

