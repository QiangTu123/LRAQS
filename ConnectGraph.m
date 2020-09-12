function [  ] = ConnectGraph( SystemModel, NeighborList, Anchor )

[m,n] = size(SystemModel);
plot(SystemModel([1:Anchor],1),SystemModel([1:Anchor],2),'r*');hold on;
plot(SystemModel([Anchor+1:m],1),SystemModel([Anchor+1:m],2),'ro');hold on;
