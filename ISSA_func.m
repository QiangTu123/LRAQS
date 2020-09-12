function [FoodPosition,FoodFitness] = ISSA_func(Beaconlocation,Distance,h,Population,Maxgen,ub,lb,Dim)

xub = ub * ones(1,Dim);
xlb = lb * ones(1,Dim);
PopulationSize = Population;

% Initialize
SalpPositions = ub .* rand(PopulationSize,Dim);
FoodFitness = inf;
for i = 1:PopulationSize
    SalpFitness(i) = Evaluation(Beaconlocation, Distance, h, SalpPositions(i,:));
end
[sorted_salps_fitness, sorted_index] = sort(SalpFitness);

Sorted_salps = zeros(PopulationSize,Dim);
for newindex = 1:PopulationSize
    Sorted_salps(newindex,:) = SalpPositions(sorted_index(newindex),:);
end
FoodPosition = Sorted_salps(1,:);
FoodFitness = sorted_salps_fitness(1);

% Main loop
loop = 2;
while loop < Maxgen+1
    c1 = 2 * exp(-(4*1/Maxgen)^2);
    mbest = mean(sum(SalpPositions(Population/2,:),1),1); % 前 N/2 个leader salps 的平均位置
    for i = 1:PopulationSize
        if i <= PopulationSize/2
            for j = 1:Dim
                c2 = rand;
                c3 = rand;
                if c3 < 0.5
                    SalpPositions(i,j) = FoodPosition(j) + c1 * ((ub-lb) * c2 + lb);
                else
                    SalpPositions(i,j) = FoodPosition(j) - c1 * ((ub-lb) * c2 + lb);
                end
            end
        else
            cn = rand;
            if cn<0.5
                point1 = SalpPositions(i-1,:);
                point2 = SalpPositions(i,:);
                SalpPositions(i,:) = (point1 + point2)/2;
            elseif cn<0.75
                b = rand;
                u = rand;
                Pa = b.*FoodPosition;
                L = c1.*abs(FoodPosition-mbest);
                SalpPositions(i,:) = Pa+(L/8)*log(exp(1/u));
            else
                b = rand;
                u = rand;
                Pa = b.*FoodPosition;
                L = c1.*abs(FoodPosition-mbest);
                SalpPositions(i,:) = Pa-(L/8)*log(exp(1/u));
            end
            
        end
    end
    for i = 1:PopulationSize
        Tp=SalpPositions(i,:)>xub;Tm=SalpPositions(i,:)<xlb;SalpPositions(i,:)=(SalpPositions(i,:).*(~(Tp+Tm)))+xub.*Tp+xlb.*Tm;
        SalpFitness(i) = Evaluation(Beaconlocation, Distance, h, SalpPositions(i,:));
        if SalpFitness(i) < FoodFitness
            FoodPosition = SalpPositions(i,:);
            FoodFitness = SalpFitness(i);
        end
    end
    loop = loop + 1;
end





