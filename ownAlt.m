function [T, a, P, rho ] = ownAlt(altitude)

%Note table now only goes up to 15k feet
altData = importdata ('altitudeTable.csv');
altTable = altData.data;

%alt  sigma  delta  theta  temp  press    dens     a    visc  k.visc  ratio

%for elements = (2:size(kGraph1))
 %   graphX = kGraph1(elements, 1);
  %  if k_XAxis < graphX && xCoordChecker == 0
   %    xCoordChecker = 1;
%  end




altChecker = 0;
for elements = (1:size (altTable))
    xOne = altTable (elements, 1);
    if altChecker == 0
        if xOne >= altitude
            altChecker = 1;
            xZero = altTable (elements-1,1);
            linInterp = (xOne - altitude) / (xOne -xZero);
            T = altTable (elements-1, 5) + (altTable (elements,5) - altTable (elements - 1,5)) * linInterp;
            T = T + 459.67; % Unit conversion from degrees rankine
            %T in deg F
            a = altTable (elements-1, 8) + (altTable (elements, 8) - altTable (elements-1, 8)) * linInterp;
            %a in ft/s
            P = altTable (elements-1, 6) + (altTable (elements, 6) - altTable (elements-1, 6)) * linInterp;
            %P in lb/ft^2
            rho = altTable (elements-1, 7) + (altTable (elements, 7) - altTable (elements-1, 7)) * linInterp;
            rho = rho * 55.93;
            %rho in lbm/m^3
            altChecker = 1;
        end
    end
end

        
        