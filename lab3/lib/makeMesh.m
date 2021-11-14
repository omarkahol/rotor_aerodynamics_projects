function helicopter = makeMesh(helicopter,N)
%Mesh the helicopter blade

helicopter.dR = (helicopter.Radius-helicopter.CutOut)/N;
space = 1:1:N;
helicopter.Mesh = helicopter.CutOut + helicopter.dR.*(space - 0.5);

end

