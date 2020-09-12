function [ NodeCoordinate ] = SystemModel_S( Square, NodeDensity, Hole )
% generate S-shaped networks
NodeNumber = round(NodeDensity * Square^2);
NodeCoordinate = zeros(NodeNumber,2);
for i = 1:NodeNumber
    TempX = rand*Square;
    TempY = rand*Square;
    
    while((TempX<50 && TempY>20 && TempY<40) || (TempX>50 && TempY>60 && TempY<80))
        TempX =rand*Square;
        TempY =rand*Square;
    end
    NodeCoordinate(i,1)=TempX;
    NodeCoordinate(i,2)=TempY;
end

