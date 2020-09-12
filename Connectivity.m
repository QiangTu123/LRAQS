function [ Flag ] = Connectivity( NeighborList )

Flag = 1;
m = size(NeighborList,1);
ConMatrix = NeighborList;

for i = 2:m
    NeighborList = NeighborList + ConMatrix^i;
end

for i = 1:m
    for j = 1:m
        if (NeighborList(i,j) == 0)
            Flag = 0;break;
        end
    end
end

