function [helicopter] = buildHelicopter(R,R0,C,NB,NP,GAMMA)
%BUILDHELICOPTER Summary of this function goes here
%   Detailed explanation goes here

helicopter = struct('Radius',R,'CutOut',R0,'Chord',C,'NumberOfBlades',NB);
helicopter.Solidity = NB*C/(pi*R);
helicopter.Area = pi*R^2;
helicopter.dR = (helicopter.Radius-helicopter.CutOut)/NP;
space = 1:1:NP;
helicopter.Mesh = helicopter.CutOut + helicopter.dR.*(space - 0.5);
helicopter.LockeNumber = GAMMA;

end

