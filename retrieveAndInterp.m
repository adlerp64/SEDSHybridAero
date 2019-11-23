function [val] = retrieveAndInterp (t,thrustTable,curveNum)
%function [thrust, propMass, CG, m_dot] = thrustTwo(t)

%altData = importdata ('altitudeTable.csv');
%altTable = altData.data;
%thrust table goes from 0 up in increments of 0.01 starting from cell
%2 although .data gets rid of that

%disp ("t");
%disp (t);
%tTest = t / 0.01;
%tTest = round(tTest);
%thrust2 = thrustTable (tTest+1, 6);


tTest = t / 0.01;
cellNum = (tTest);
%disp ("firstCellNum");
%disp (cellNum);
cellNum = cellNum + 0.5;
cellNumRound = round (cellNum);


%Linear Interpolation
%tSmall = thrustTable (cellNumRound, 1);
%tBig = thrustTable (cellNumRound+1,1);
disp (cellNumRound);
val = thrustTable (cellNumRound, curveNum);
%val =  (thrustTable(cellNumRound+1, curveNum) - thrustTable (cellNumRound, curveNum)) / 0.01 + thrustTable (cellNumRound, curveNum);




