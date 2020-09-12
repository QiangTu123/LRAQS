function [ NodeCoordinate ] = SystemModel_O( Square, NodeDensity, Hole )
% generate S-shaped networks
NodeNumber = round(NodeDensity * Square^2);
NodeCoordinate = zeros(NodeNumber,2);
for i = 1:NodeNumber
    TempX = rand*Square;
    TempY = rand*Square;
    
    while((TempX>(Square/2-Hole/2) && TempX<(Square/2+Hole/2)) && (TempY>(Square/2-Hole/2) && TempY<(Square/2+Hole/2)))
        TempX =rand*Square;
        TempY =rand*Square;
    end
    NodeCoordinate(i,1)=TempX;
    NodeCoordinate(i,2)=TempY;
end


