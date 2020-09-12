%=========================================================================%
% Title:Range-free localization using reliable anchor pair selection and quantum-behaved
%       salp swarm algorithm for anisotropic wirless sensor networks (LRAQS£©
% Developed in MATLAB R2014a
%=========================================================================%
close all;
clc;clear;

%------------------------------Parameter Initialization-----------------------------%
Square = 100;      % Side length of square area
CommRange = 20;    % Communication radii of unknown node
AnCommRange = 20;  % Communication radii of anchor node
Hole = 40;         % Side length of hole 
DoI = 0.05;        % Degree of Irregularity
Trial = 1;       % Number of trials
AvgErr = zeros(1,Trial);
AvgTest2 = zeros(1,Trial);
AnchorL = 10:5:50;                % Number of anchor nodes
% Density = 0.006:0.002:0.03;     % Node density
% Density = 0.006;
ERROR = zeros(1,length(AnchorL)); % Mean localization error
MDE = zeros(1,length(AnchorL));   % Average distance error
MaxHop = 4; % Threshold of maximum hop counts 

for De = 1:length(AnchorL)  % Different number of anchor nodes
    NodeDensity = 0.01;  % Current node density
    Anchor = AnchorL(De);  % Current number of anchor nodes
    for Tr = 1:Trial
        close all;
        disp('The Number of Anchor = ');disp(Anchor);
        disp('Trial = ');disp(Tr);
        %%%% SystemModel_O %%%%%%%
        NodeCoordinate = SystemModel_O(Square,NodeDensity,Hole); 
        NeighborList = Neighbor(NodeCoordinate,CommRange,AnCommRange,Anchor,DoI); 

        while(~Connectivity(NeighborList))
            NodeCoordinate = SystemModel_O(Square,NodeDensity,Hole);
            NeighborList = Neighbor(NodeCoordinate,CommRange,AnCommRange,Anchor,DoI);
        end
       %% ~~~~~~~~~~~~ Network Parameter Initialization ~~~~~~~~~~~%%
        NodeAmount = size(NeighborList,1); % Number of sensor nodes
        UnAmount = NodeAmount - Anchor;    % Number of unknown nodes
        BeaconLocation = NodeCoordinate(1:Anchor,:); % Location of anchor nodes
        UnLocation = NodeCoordinate(Anchor+1:NodeAmount,:); % Real location of unknown nodes
        TP = zeros(UnAmount,2); % Estimated location of unknown nodes
        NeighborArray = NeighborList; 
        Dall = zeros(NodeAmount,NodeAmount); 
        
        ConnectGraph(NodeCoordinate,NeighborList,Anchor);
        
        % Initialize the matrix of hop count 
        for i = 1:size(NeighborList,1)
            for j = 1:size(NeighborList,2)
                Dall(i,j) = sqrt((NodeCoordinate(i,1)-NodeCoordinate(j,1))^2 + (NodeCoordinate(i,2)-NodeCoordinate(j,2))^2);
                if NeighborList(i,j) == 0
                    NeighborList(i,j) = inf;
                end
            end
        end
        NeighborList(logical(eye(size(NeighborList,1)))) = 0;
        
        for k=1:size(NeighborList,1)
            for i=1:size(NeighborList,1)
                for j=1:size(NeighborList,1)
                    if NeighborList(i,k)+NeighborList(k,j)<NeighborList(i,j) %min(h(i,j),h(i,k)+h(k,j))
                        NeighborList(i,j)=NeighborList(i,k)+NeighborList(k,j);
                    end
                end
            end
        end
        
        % Compute the expected hop progress of network
        an = pi/4; % [0,pi/2]
        % Equation.(18)
        hs_e = 2 * NodeDensity * sin(an) * quadgk(@(L)((L.^2) .* exp(-NodeDensity * (an) * CommRange^2 + NodeDensity * (an) * L.^2)),0,CommRange);
        hopBtoU = NeighborList(1:Anchor,(Anchor+1:NodeAmount)); 
        Distance = hs_e .* hopBtoU;
        
       %% Use QSSA to calculate the location of unknown nodes
        % Parameter Initialization of QSSA
        PopulationSize = 20;  % Population size
        MaxGen = 200;  % Maximum of iterations
        PreLocation = zeros(1,2); % Estimated location of one unknown node
        
        %% Correct the estimated distance between anchor node and unknown node based on RAPS
        NodeQueList = NeighborList(Anchor+1:NodeAmount,1:Anchor);
        ReliaAnSet = zeros(Anchor); % Construct the set of reliability of anchor pair
        NodePosition = zeros(size(NodeQueList,1),2); %  Location of reliable anchor nodes
        D1 = Dall(1:Anchor,1:Anchor); % Distance matrix of anchors
        
        for k = 1:UnAmount
            % Compute the reliability of each anchor pairs
            for i = 1:Anchor
                for j = 1:Anchor
                    ReliaAnSet(i,j) = D1(i,j)/(hopBtoU(i,k)+hopBtoU(j,k)); % Eq.(4)
                end
            end
            %
            DisList = zeros(1, Anchor); %  Distance between unkonwn node and anchor node 
            for i = 1:Anchor
                In = 0; 
                while (DisList(i) == 0 && In< Anchor/2)
                    j = find(ReliaAnSet(i,:)==max(ReliaAnSet(i,:))); 
                    In = In+1;
                    if i ~= j
                        TempRelia = ReliaAnSet(i,j);
                        ReliaAnSet(i,j) = 0; 
                        D = D1(i,j);
                        % Use eq.(14) to compute the distance between unkonwn node k and anchor j if the pair is super pair, 
                        if ((((NodeQueList(k,i)*CommRange)^2 + D^2 - (NodeQueList(k,j)*CommRange)^2)/(2*NodeQueList(k,i)*CommRange*D)<=1 &&... %% Super pair
                                ((NodeQueList(k,i)*CommRange)^2 + D^2 - (NodeQueList(k,j)*CommRange)^2)/(2*NodeQueList(k,i)*CommRange*D)>= -1)&&((NodeQueList(k,i)*CommRange)<D &&( NodeQueList(k,j)*CommRange)<D))
                            DisList(i) = NodeQueList(k,i)*TempRelia*real(Integral1(NodeQueList(k,i), NodeQueList(k,j), CommRange, D));
                            Distance(i,k) = DisList(i);
                        end
                    end
                end
            end
      
            TempNodeQueList = NodeQueList(k,:);
            TempBeaconLocation = BeaconLocation;
            TemphopBtoU = hopBtoU;
            TempDistance = Distance;
            
           % Utilize maxhop to limit the number of transmission between nodes
            CastList = find(TempNodeQueList>MaxHop); 
            TempBeaconLocation(CastList,:) = [;];
            TemphopBtoU(CastList,:) = [;];
            TempDistance(CastList,:) = [;];
            
            % Use QSSA to replace the maximum likelihood estimated method
            [PreLocation, fitnessvalue] = ISSA_func(TempBeaconLocation, TempDistance(:,k),TemphopBtoU(:,k),PopulationSize,MaxGen,Square,0,2);
            TP(k,1) = PreLocation(1,1);
            TP(k,2) = PreLocation(1,2);
        end
        
       %% Metrics: ERROR and MDE
        ER = 0; % Mean localization error
        for Node = 1:UnAmount
            ER = ER + sqrt((UnLocation(Node,1)-TP(Node,1))^2+(UnLocation(Node,2)-TP(Node,2))^2);
        end
        
        DR = 0; % Mean distance error
        DPredict = zeros(Anchor,UnAmount); 
        DTure = Dall(1:Anchor,Anchor+1:NodeAmount); 
        for i = 1:Anchor
            for j = 1:UnAmount
                DPredict(i,j) = sqrt((BeaconLocation(i,1)-TP(j,1))^2+(BeaconLocation(i,2)-TP(j,2))^2);
            end
        end
        DR = (sum((sum((abs(DPredict-DTure)),2)),1))/(Anchor*UnAmount); 
        MDE(De) = MDE(De) + DR;
        
        disp('Localization error ='); disp(ER/UnAmount);  % Localization error
        disp('Distance error'); disp(DR) % Distance error
    end
    
    ERROR(De) = ERROR(De)/(Trial);
    MDE(De) = MDE(De)/Trial;
    disp(ERROR);
    disp(MDE);
end

