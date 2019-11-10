%function [suc, transYieldLoad] = transYieldFailure (Fy, D, t, W, a, alum, ultYieldStress, margin)
%Damping force : dampForce
%	Input: []
%	Output: (dampForce)
function [dampTorq, dampCof] = dampTorque (mDot, dCOM, omega)

%Distances from nose cone base to elements
inToM = 0.0254;
lbToKg = 0.4536;

dNozzle = 10 * 12 * inToM;

dCOMOx = 6.5 * 12 * inToM;
dCOMProp = 9.5 * 12 * inToM;

mOx = 21;
mProp = 3;

dCOM_fuel = (dCOMOx * mOx + dCOMProp * mProp) / (mOx + mProp);

dCOM_Exit = dNozzle - dCOM;
dCOM_fuel = abs(dCOM - dCOM_fuel);

dampCof = mDot * ( dCOM_Exit^2 - dCOM_fuel^2);
if (dCOM_fuel == 0)
    dampCof = 0;
end


dampTorqMetric = dampCof * omega;
dampTorq = dampTorqMetric * 1/inToM * 1/lbToKg;
dampCof = dampCof * 1/inToM * 1/lbToKg;
