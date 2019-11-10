

%function [drag] = dragForce((phi-alpha),rho,v)
%function [normalForceExp, cpLoc] = normalForce(alpha, rho, v)

%function [torque] = torqueFromForce (angleForce, magnitudeForce, angleLength, magnitudeLength)
%momLength = (COPDif^2 + xDev^2)^0.5;
%momAngle = alpha - atan(xDev/COPDif); %-aTan because counterclockwise is positive

rho = 0.0408;
xDev = 0.02/12;
cog = 5.0826;

% v: 400-750 ft/s
% phi-alpha: 0-5: increments of 0.1


comparisons = [1,2,3];
for v = 400:50:750
    newRow = [v,v,v];
    comparisons = [comparisons;newRow];
    for angle=0.1:0.1:5
        angleRad = angle * 0.01745;
        drag = 0.6 * 3.14159 * (3/12)^2 * rho * v^2;
        [fN, cpLoc] = normalForce (angleRad, rho, v);
        cpDif = abs (cog - cpLoc);
        momLength = (cpDif^2 + xDev^2)^0.5;
        momAngle = angleRad - atan (xDev/ cpDif);
        
        dragTorque = torqueFromForce(3.14159/2 + angleRad, drag, momAngle, momLength);
        normalTorque = torqueFromForce (3.14159/2, fN, 0, cpDif);
        newRow = [angle, dragTorque, normalTorque];
        comparisons = [comparisons;newRow];
    end
end
csvwrite ('dragTNormalTComp.csv', comparisons);
        
        
        

