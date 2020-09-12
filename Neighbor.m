function [ NeighborList ] = Neighbor( SystemModel, CommRange, AnCommRange, Anchor, DoI )

[m,n] = size(SystemModel);
NeighborList = zeros(m);

for i = 1:Anchor
    for j = i:m
        if (sqrt((SystemModel(i,1)-SystemModel(j,1))^2 + (SystemModel(i,2)-SystemModel(j,2))^2)/AnCommRange <= (1-DoI))
            NeighborList(i,j) = 1;
            NeighborList(j,i) = 1;
        elseif (sqrt((SystemModel(i,1)-SystemModel(j,1))^2 + (SystemModel(i,2)-SystemModel(j,2))^2)/AnCommRange >= (1+DoI))
            NeighborList(j,i) = 0;
            NeighborList(i,j) = 0;
        else
            P = (sqrt((SystemModel(i,1)-SystemModel(j,1))^2 + (SystemModel(i,2)-SystemModel(j,2))^2)-AnCommRange*(1-DoI))/(2*AnCommRange*DoI);
            NeighborList(j,i) = sign(fix(rand/P));
            NeighborList(i,j) = NeighborList(j,i);
        end
    end
end

for i = Anchor + 1:m
    for j = i:m
        if (sqrt((SystemModel(i,1)-SystemModel(j,1))^2 + (SystemModel(i,2)-SystemModel(j,2))^2)/CommRange <= (1-DoI))
            NeighborList(i,j) = 1;
            NeighborList(j,i) = 1;
        elseif (sqrt((SystemModel(i,1)-SystemModel(j,1))^2 + (SystemModel(i,2)-SystemModel(j,2))^2)/CommRange >= (1+DoI))
            NeighborList(j,i) = 0;
            NeighborList(i,j) = 0;
        else
            P = (sqrt((SystemModel(i,1)-SystemModel(j,1))^2 + (SystemModel(i,2)-SystemModel(j,2))^2)-CommRange*(1-DoI))/(2*CommRange*DoI);
            NeighborList(j,i) = sign(fix(rand/P));
            NeighborList(i,j) = NeighborList(j,i);
        end
    end
end

for i = 1:m
    NeighborList(i,i) = 0;
end

