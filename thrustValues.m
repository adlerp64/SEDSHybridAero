%function [thrust, propMass, CG, m_dot] = thrustTwo(t)
function [thrust, propMass, CG, m_dot] = thrustValues (t, thrustTable)


m_prop = 23.882;
m_dry = 95;                % [lbm] dry mass of rocket 
x_dry = 5.082569431;        % [lbm] dry CG location
x_tank = 6.5;               % [ft] location of middle of oxidizer tank
x_fg = 9.5;                 % [ft] location of middle of fuel grain
x_prop = 0.5*(x_tank + x_fg);

%function [val] = retrieveAndInterp (t,thrustTable,curveNum)
if t < 9.50
    thrust = retrieveAndInterp (t, thrustTable, 6);
    m_dot = retrieveAndInterp (t,thrustTable, 7);
    propMass = m_prop - (retrieveAndInterp(t,thrustTable,13));
    
else
    thrust = 0;
    m_dot = 0;
    propMass = 0;
end


%m_prop = 27; % [lbm] propellant mass % Updated value is 24.


CG = (m_dry*x_dry + propMass*x_prop) / (m_dry + m_prop);


