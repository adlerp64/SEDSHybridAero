%function [k] = sortGraph (k_XAxis, k_YAxis, kGraph1, secondAxis)
function [dDampT] = dragDamp (omega, rho, com, l)

%F = 1/2 * Cd*rho*v^2*L*D

%Better equation for this based on velocity?
Cd = 0.2;

%These were from integrals of drag equation from tip to Cd and tail to Cd
dDampUpper = 1/6*rho*omega*((l-com))^4;
dDampLower = 1/6*rho*omega*(com)^4;

%-1 is to see it go opposite direction of omega
dDampT = -1 * Cd* (dDampLower + dDampUpper);

