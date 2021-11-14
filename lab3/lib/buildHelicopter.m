function [helicopter] = buildHelicopter(R,R0,c,Nb)
%BUILDHELICOPTER Summary of this function goes here
%   Detailed explanation goes here

helicopter = struct('Radius',R,'CutOut',R0,'Chord',c,'NumberOfBlades',Nb);
helicopter.Solidity = Nb*c/(pi*R);
helicopter.Area = pi*R^2;
helicopter.Mesh = [];
helicopter.dR = 0;

end

