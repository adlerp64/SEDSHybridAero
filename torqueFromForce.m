function [torque] = torqueFromForce (angleForce, magnitudeForce, angleLength, magnitudeLength)

%AngleForce is the angle couterclockwise from the direction of gravity
%(AngleForce is starting from 90 degrees of traditional unit circle)

%^Same with angleLength

fy = cos(angleForce) * magnitudeForce;
fx = -1 * sin(angleForce) * magnitudeForce;

Ly = cos(angleLength) * magnitudeLength;
Lx = -1 * sin(angleLength) * magnitudeLength;

torque = fy*Lx - fx*Ly;