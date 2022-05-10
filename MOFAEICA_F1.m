%% MOFAEICA for F1 (MOP1)
% clc;
clear;
close all;
%MaxRun = 30;
MaxRun = 1; % For the capsule
%MaxFEs = 40000;
MaxFEs = 2100;
FEsRuns = zeros(MaxRun, 1);
IGDRuns= zeros(MaxRun, 1);
NHVRuns= zeros(MaxRun, 1);
FEs = 0;
for nRun=1:MaxRun
disp('MOFAEICA on F1(MOP1)...');
disp(['To reduce compute time, MaxFEs was set to ' num2str(MaxFEs)]);
disp(['Independant Run : ' num2str(nRun)]);
%% Problem Definition
CostFunction = @F1; % MOP1
nVar = 30;
down=[0]; %#ok<*NBRAK>
up=[+1];
VarSize = [1 nVar];
% Initial settings for MOFAEICA
MaxIt = 9472; % Stopping criterion: FEs at each iteration*MaxIt == 40,000 for F1, so MaxIt can be any number larger than 40,000/FEs at each iteration
Alpha = 10; % This Alpha is used for the change range in the velocity component of each solution
nPop = 250;
nEmp = 50;
nCol = nPop - nEmp;

C=0.3;
% Probability of applying Velocity Divergence
PVD = 1;

% Probability of applying Velocity Adaptation
PVA = 1;

% Probability of applying Information Sharing between empires
PLS = 1;

for j = 1:nVar
  GlobalBest.Position(j) = rand.*(up-down) + down;
end
GlobalBest.Velocity = zeros(VarSize);
GlobalBest.Cost = Inf;
GlobalBest.Power = 0; % The worst value

%% Initialization
empty_country.Position=[];
empty_country.Velocity=[];
empty_country.Cost=[];
empty_country.Power=[];
empty_country.Rank=[];
empty_country.DominationSet=[];
empty_country.DominatedCount=[];
empty_country.PS=[];

empty_country.Best.Position=[];
empty_country.Best.Velocity=[];
empty_country.Best.Cost=[];
empty_country.Best.Power=[];
empty_country.Best.Rank=[];
empty_country.Best.DominationSet=[];
empty_country.Best.DominatedCount=[];
empty_country.Best.PS=[];

empty_Archive.Position = [];
empty_Archive.Velocity = [];
empty_Archive.Cost = [];
empty_Archive.Power = [];
empty_Archive.Rank = [];
empty_Archive.DominationSet = [];
empty_Archive.DominatedCount = [];
empty_Archive.PS = [];

empty_Archive.Best.Position = [];
empty_Archive.Best.Velocity = [];
empty_Archive.Best.Cost = [];
empty_Archive.Best.Power = [];
empty_Archive.Best.Rank = [];
empty_Archive.Best.DominationSet = [];
empty_Archive.Best.DominatedCount = [];
empty_Archive.Best.PS = [];

country=repmat(empty_country,nPop,1);

BestPower = zeros(MaxIt,1);
BestCost = zeros(MaxIt,1);
GlobalBestCost = zeros(1, 1);
GlobalBestPosition = zeros(1, nVar);
convergence1 = zeros(MaxIt, 1);

  for i=1:nPop
      for j = 1:nVar
          country(i).Position(j) = rand.*(up-down) + down;
      end      
      country(i).Velocity = zeros(VarSize);
      country(i).Cost = CostFunction(country(i).Position);      
      country(i).Power = 0;

      country(i).Best.Position = country(i).Position;
      country(i).Best.Velocity = country(i).Velocity;                 
      country(i).Best.Cost = country(i).Cost;      
      country(i).Best.Power = country(i).Power;            
      country(i).Best.Rank = country(i).Rank;                                                                                                                   
  end
% Fast Non-Dominated Sorting
%fprintf('Fast Non-dominated sorting, Penalized Sigma Index for Fitness (Power) estimation, and Sub-population forming ... \n');
[country, F] = NonDominatedSorting(country);
% Calculate Density Estimation in objective and decision spaces
country = PS(country,F);
% Calculate Fitness (Power)
country = PowerEstimator(country);

% Update Best Experience and Global Best Solution
for i = 1:nPop
      country(i).Best.Position = country(i).Position;
      country(i).Best.Velocity = country(i).Velocity;           
      country(i).Best.Cost = country(i).Cost;      
      country(i).Best.Power = country(i).Power;            
      country(i).Best.Rank = country(i).Rank;
      country(i).Best.DominationSet = country(i).DominationSet;
      country(i).Best.DominatedCount = country(i).DominatedCount;
      country(i).Best.PS = country(i).PS;
      
      if country(i).Best.Power > GlobalBest.Power
              GlobalBest = country(i).Best;
      end 
end

%%  Form imperialists and colonies    
% Sort countries
[~,index] = sort([country.Power],'descend');
country=country(index);
% Assign of Colonies and Imperialists
imp=country(1:nEmp);
col=country(nEmp+1:nPop);

empty_empire.Imp=[];
empty_empire.Col=repmat(empty_country,0,1);
empty_empire.nCol=0;
emp=repmat(empty_empire,nEmp,1);

% Assign Imperialists
    for k=1:nEmp
        emp(k).Imp=imp(k);
    end
        
% Assign Colonies
Np = [imp.Power]./max([imp.Power]);
SNp = sum(Np);
Npp = Np./SNp;

for Num=1:nCol
     k = RouletteWheelSelection(Npp);
     check = isempty(emp(k).Col);
     if check == 1
         emp(k).Col = col(Num);
     else
        emp(k).Col = [emp(k).Col; col(Num)];
     end
     emp(k).nCol = emp(k).nCol + 1;                 
end

index = find([emp.nCol]==0);
if index>=1
    for tedademp=1:numel(index)
        [quan, ind]=max([emp.nCol]);
        randcol = randi(quan);
        emp(index(tedademp)).Col = emp(ind).Col(randcol);
        emp(index(tedademp)).nCol = 1;
        emp(ind).Col(randcol) = [];
        emp(ind).nCol = emp(ind).nCol - 1;             
    end  
end
Powers = zeros (MaxIt, 1);
Archive = repmat(empty_Archive,1,1);

 %% Main Loop   
for it=1:MaxIt
fprintf('Iteration time is %d \n', it)
% FIS1 is used by Global Learning, Universal global best diversity,
% Differential evolutionary-based local search to adapt the parameters at
% each iteration time

% FIS2 is used for dynamic selection of operators at each time window
if mod(it, MaxFEs/10) == 0 % tw = MaxFEs/10;
%% Fuzzy Adaptive Operator Selection
    BestPower_tw = BestPower(it-10+1:it-1);
    Delta_BestPower_tw =  max(BestPower_tw) - min(BestPower_tw);
    Stagnation = 1 - (Delta_BestPower_tw / max(BestPower_tw));
    Stagnation = min(Stagnation, 1);
    Stagnation = max(Stagnation, 0);
    UFuzzy = [Stagnation; PVA; PVD; PLS];
    FISMAT = readfis('FAOS.fis');
    Y = evalfis(FISMAT, UFuzzy);
    PVA = Y(1,1);
    PVD = Y(1,2);
    PLS = Y(1,3);
end

%% Global learning–based velocity adaptation
if rand <= PVA
for k = 1:numel(emp)
    for col = 1:numel(emp(k).Col)
        NP = (abs(emp(k).Imp.Power - emp(k).Col(col).Power))/GlobalBest.Power;
        NP = min(NP,1);
        NP = max(NP,0);
        NFEs = FEs/MaxFEs;            
        UFuzzy = [NP; NFEs];
        FISMAT = readfis('FAGLVA.fis');
        Y = evalfis(FISMAT, UFuzzy);
        % Social learning parameters {w, c2, Beta}
        % Cognitive learning parameters {c1}
        w = Y(1,1);
        c1 = Y(1,2);
        Beta = Y(1,3);
        c2 = Y(1,4);

        tempCol1 = emp(k).Col(col);
        tempImp1 = emp(k).Imp;                                   
        emp(k).Col(col).Velocity = (w.*emp(k).Col(col).Velocity)+(Beta.*rand(VarSize)).*(emp(k).Imp.Position-emp(k).Col(col).Position) + (c1.*rand(VarSize)).*(emp(k).Col(col).Best.Position-emp(k).Col(col).Position) + (c2.*rand(VarSize)).*(GlobalBest.Position-emp(k).Col(col).Position);

        [VelMin, VelMax] = VelLimit(GlobalBest.Position, emp(k).Col(col).Position, it, up, down, Alpha);
        emp(k).Col(col).Velocity = min(max(emp(k).Col(col).Velocity,VelMin),VelMax);
        [VelMin, VelMax] = VelLimit(GlobalBest.Position, emp(k).Imp.Position, it, up, down, Alpha);
        emp(k).Imp.Velocity = min(max(emp(k).Imp.Velocity,VelMin),VelMax);                   
        emp(k).Col(col).Position = emp(k).Col(col).Position + emp(k).Col(col).Velocity;
        emp(k).Imp.Position = emp(k).Imp.Position + emp(k).Imp.Velocity;
        for flg=1:nVar
            if emp(k).Col(col).Position(flg) < down || emp(k).Col(col).Position(flg) > up 
                emp(k).Col(col).Velocity(flg) = -emp(k).Col(col).Velocity(flg);              
            end
        end
        for flg=1:nVar
            if emp(k).Imp.Position(flg) < down || emp(k).Imp.Position(flg) > up %#ok<*BDSCI>
                emp(k).Imp.Velocity(flg) = -emp(k).Imp.Velocity(flg);              
            end
        end
        for flg=1:nVar
                emp(k).Imp.Position(flg) = min(up,max(down,emp(k).Imp.Position(flg))); % Bound the new location
        end
        for flg=1:nVar
                emp(k).Col(col).Position(flg) = min(up,max(down,emp(k).Col(col).Position(flg))); % Bound the new location
        end

        emp(k).Col(col).Cost = CostFunction(emp(k).Col(col).Position);
        FEs = FEs + 1;            

%         if emp(k).Col(col).Power >= emp(k).Col(col).Best.Power
%             emp(k).Col(col).Best.Position = emp(k).Col(col).Position;
%             emp(k).Col(col).Best.Velocity = emp(k).Col(col).Velocity;
%             emp(k).Col(col).Best.Cost = emp(k).Col(col).Cost;
%             emp(k).Col(col).Best.Power = emp(k).Col(col).Power;
%         end

%         if emp(k).Col(col).Power >= GlobalBest.Power
%             temp = GlobalBest;
%             GlobalBest.Position = emp(k).Col(col).Position;
%             GlobalBest.Velocity = emp(k).Col(col).Velocity;
%             GlobalBest.Cost = emp(k).Col(col).Cost;
%             GlobalBest.Power = emp(k).Col(col).Power;
%             emp(k).Col(col).Position = temp.Position;
%             emp(k).Col(col).Velocity = temp.Velocity;
%             emp(k).Col(col).Cost = temp.Cost;
%             emp(k).Col(col).Power = temp.Power;
%         end

        % Exchange a colony with its local imperialist if the colony has more power than its local imperialist
%         if emp(k).Col(col).Power > emp(k).Imp.Power
%             [emp(k).Imp, emp(k).Col(col)] = deal(emp(k).Col(col), emp(k).Imp);                
%         end

        emp(k).Imp.Cost = CostFunction(emp(k).Imp.Position);
        FEs = FEs + 1;

%         if emp(k).Imp.Power >= emp(k).Imp.Best.Power
%             emp(k).Imp.Best.Position = emp(k).Imp.Position;
%             emp(k).Imp.Best.Velocity = emp(k).Imp.Velocity;
%             emp(k).Imp.Best.Cost = emp(k).Imp.Cost;
%             emp(k).Imp.Best.Power = emp(k).Imp.Power;
%         end    

%         if emp(k).Imp.Power >= GlobalBest.Power
%             temp = GlobalBest;
%             GlobalBest.Position = emp(k).Imp.Position;
%             GlobalBest.Velocity = emp(k).Imp.Velocity;
%             GlobalBest.Cost = emp(k).Imp.Cost;
%             GlobalBest.Power = emp(k).Imp.Power;
%             emp(k).Imp.Position = temp.Position;
%             emp(k).Imp.Velocity = temp.Velocity;
%             emp(k).Imp.Cost = temp.Cost;
%             emp(k).Imp.Power = temp.Power;
%         end        
    end
end
end

%% Universal global best diversity
if rand <= PVD
for k = 1:numel(emp)
    for col = 1:numel(emp(k).Col)
        NP = abs(emp(k).Imp.Power - emp(k).Col(col).Power)/GlobalBest.Power;
        NP = min(NP,1);
        NP = max(NP,0);
        NFEs = FEs/MaxFEs;             
        UFuzzy = [NP; NFEs];
        FISMAT = readfis('FAUDVD.fis');
        Y = evalfis(FISMAT, UFuzzy);        
        Pdiv = Y(1,1);
        Np = [imp.Power]./max([imp.Power]);
        SNp = sum(Np);     
        if SNp == 0
            Npp = Np;
        else
        Npp = Np./SNp;
        end             

        tempCol = emp(k).Col(col);
        for d = 1:nVar
            randk = RouletteWheelSelection(Npp);
            while numel(emp(randk).Col) == 0
                randk = RouletteWheelSelection(Npp);
            end
            Npcol = [emp(randk).Col.Power]./max([emp(randk).Col.Power]);
            SNpcol = sum(Npcol);
            if SNpcol == 0
                Nppcol = Npcol;
            else
            Nppcol = Npcol./SNpcol;            
            end
            status = isnan(Nppcol);
            if (Nppcol ~= 0) & (status ~= 1) %#ok<*AND2>
                randcol = RouletteWheelSelection(Nppcol);
            else
                randcol = randi(numel(emp(randk).Col));
            end
            emp(k).Col(col).Velocity(d) = emp(k).Col(col).Velocity(d) + rand.*(emp(randk).Col(randcol).Best.Position(d) - emp(k).Col(col).Position(d));
        end

        tempGL = GlobalBest;
        for d = 1:nVar
            randk = RouletteWheelSelection(Npp);
            if rand < PVD
                GlobalBest.Velocity(d) = GlobalBest.Velocity(d) + rand.*(emp(randk).Imp.Best.Position(d) - GlobalBest.Position(d));
            end
        end

        % Constraints are applied to the velocities
        [VelMin, VelMax] = VelLimit(GlobalBest.Position, emp(k).Col(col).Position, it, up, down, Alpha);
        emp(k).Col(col).Velocity = min(max(emp(k).Col(col).Velocity,VelMin),VelMax);           
        [VelMin, VelMax] = VelLimit(GlobalBest.Position, emp(k).Imp.Position, it, up, down, Alpha);
        emp(k).Imp.Velocity = min(max(emp(k).Imp.Velocity,VelMin),VelMax);
        VelMax = 0; VelMin = 0;           

        GlobalBest.Velocity = min(max(GlobalBest.Velocity,VelMin),VelMax);

        emp(k).Col(col).Position = emp(k).Col(col).Position + emp(k).Col(col).Velocity;
        emp(k).Imp.Position = emp(k).Imp.Position + emp(k).Imp.Velocity;
        GlobalBest.Position = GlobalBest.Position + GlobalBest.Velocity;

        for flg=1:nVar
            if emp(k).Col(col).Position(flg) < down | emp(k).Col(col).Position(flg) > up %#ok<*OR2>
                emp(k).Col(col).Velocity(flg) = -emp(k).Col(col).Velocity(flg);              
            end
        end
        for flg=1:nVar
                emp(k).Col(col).Position(flg) = min(up,max(down,emp(k).Col(col).Position(flg))); % Bound the new location
        end

        for flg=1:nVar
            if emp(k).Imp.Position(flg) < down | emp(k).Imp.Position(flg) > up
                emp(k).Imp.Velocity(flg) = -emp(k).Imp.Velocity(flg);              
            end
        end
        for flg=1:nVar
                emp(k).Imp.Position(flg) = min(up,max(down,emp(k).Imp.Position(flg))); % Bound the new location
        end

        for flg=1:nVar
            if GlobalBest.Position(flg) < down | GlobalBest.Position(flg) > up
                GlobalBest.Velocity(flg) = -GlobalBest.Velocity(flg);              
            end
        end
        for flg=1:nVar
                GlobalBest.Position(flg) = min(up,max(down,GlobalBest.Position(flg))); % Bound the new location
        end
        emp(k).Col(col).Cost = CostFunction(emp(k).Col(col).Position);
        FEs = FEs + 1;

%         % Update personal best of colonies
%         if emp(k).Col(col).Power >= emp(k).Col(col).Best.Power
%             emp(k).Col(col).Best.Position = emp(k).Col(col).Position;
%             emp(k).Col(col).Best.Velocity = emp(k).Col(col).Velocity;
%             emp(k).Col(col).Best.Cost = emp(k).Col(col).Cost;
%             emp(k).Col(col).Best.Power = emp(k).Col(col).Power;
%         end

%         % Update the Global Best imperialist using colonies
%         if emp(k).Col(col).Power >= GlobalBest.Power
%             temp = GlobalBest;
%             GlobalBest.Position = emp(k).Col(col).Position;
%             GlobalBest.Velocity = emp(k).Col(col).Velocity;
%             GlobalBest.Cost = emp(k).Col(col).Cost;
%             GlobalBest.Power = emp(k).Col(col).Power;
%             emp(k).Col(col).Position = temp.Position;
%             emp(k).Col(col).Velocity = temp.Velocity;
%             emp(k).Col(col).Cost = temp.Cost;
%             emp(k).Col(col).Power = temp.Power;
%         end

%         % Exchange a colony with its local imperialist if the colony has more power than its local imperialist
%         if emp(k).Col(col).Power > emp(k).Imp.Power
%             [emp(k).Imp, emp(k).Col(col)] = deal(emp(k).Col(col), emp(k).Imp);
%         end        

        emp(k).Imp.Cost = CostFunction(emp(k).Imp.Position);
        FEs = FEs + 1;

%        % Update personal best of the imperialist1
%         if emp(k).Imp.Power >= emp(k).Imp.Best.Power
%             emp(k).Imp.Best.Position = emp(k).Imp.Position;
%             emp(k).Imp.Best.Velocity = emp(k).Imp.Velocity;
%             emp(k).Imp.Best.Cost = emp(k).Imp.Cost;
%             emp(k).Imp.Best.Power = emp(k).Imp.Power;
%         end

%         % Update the Global Best imperialist using colonies
%         if emp(k).Imp.Power >= GlobalBest.Power
%             temp = GlobalBest;
%             GlobalBest.Position = emp(k).Imp.Position;
%             GlobalBest.Velocity = emp(k).Imp.Velocity;
%             GlobalBest.Cost = emp(k).Imp.Cost;
%             GlobalBest.Power = emp(k).Imp.Power;
%             emp(k).Imp.Position = temp.Position;
%             emp(k).Imp.Velocity = temp.Velocity;
%             emp(k).Imp.Cost = temp.Cost;
%             emp(k).Imp.Power = temp.Power;
%         end        

        GlobalBest.Cost = CostFunction(GlobalBest.Position);
        FEs = FEs + 1;            
%         % Global Best is sensitive, selection phase is required!       
%         if tempGL.Power > GlobalBest.Power
%              [GlobalBest, tempGL] = deal(tempGL, GlobalBest);                
%         end                         
    end
end   
end

%% Differential evolutionary-based local search
if rand <= PLS
    for k = 1:numel(emp)
        if numel(emp(k).Col) ~= 0
                temp = [emp(k).Col.Power];
                [~, indexwCol] = min(temp);
                wCol = emp(k).Col(indexwCol);   
                line_found = find([emp.nCol]~=0);
                rand1 = line_found(randi(numel(find([emp.nCol]~=0))));

                Colr1 = emp(rand1).Col(randi(numel(emp(rand1).Col)));
                    if Colr1.Position == wCol.Position
                        line_found = find([emp.nCol]~=0);
                        rand1 = line_found(randi(numel(find([emp.nCol]~=0))));                            
                        if numel(emp(rand1).Col) == 0
                            line_found = find([emp.nCol]~=0);
                            rand1 = line_found(randi(numel(find([emp.nCol]~=0))));
                        end
                    end
                Colr1 = emp(rand1).Col(randi(numel(emp(rand1).Col)));

                line_found = find([emp.nCol]~=0);
                rand2 = line_found(randi(numel(find([emp.nCol]~=0))));                   
                Colr2 = emp(rand2).Col(randi(numel(emp(rand2).Col)));     
                if Colr2.Position == wCol.Position
                    line_found = find([emp.nCol]~=0);
                    rand2 = line_found(randi(numel(find([emp.nCol]~=0))));

                    if numel(emp(rand2).Col) == 0
                        line_found = find([emp.nCol]~=0);
                        rand2 = line_found(randi(numel(find([emp.nCol]~=0))));                            
                    end                            
                    Colr2 = emp(rand2).Col(randi(numel(emp(rand2).Col)));                        
                end

                % The NP is calculated for the Colr3 and Imp
                NP = abs(emp(k).Imp.Power - wCol.Power)/GlobalBest.Power;
                NP = min(NP,1);
                NP = max(NP,0);
                NFEs = FEs/MaxFEs;
                UFuzzy = [NP; NFEs];
                FISMAT = readfis('FADELS.fis');
                Y = evalfis(FISMAT, UFuzzy);   
                FF1 = Y(1,1); %F1 has been previously used as a function, to avoid conflicting with that, we use FF1 as the name of the variable
                FF2 = Y(1,2);
                pCR = Y(1,3);
                MutantCol.Velocity = (rand*FF1).*(emp(k).Imp.Position - wCol.Position) + (rand*FF2).*(Colr1.Position - Colr2.Position);

                % Constraints are applied to the velocities            
                VelMax = +Alpha.*((up-down)./up); VelMin = -VelMax;                     
                MutantCol.Velocity = min(max(MutantCol.Velocity,VelMin),VelMax);

                for d=1:nVar
                    if rand <= PLS
                        TrialCol.Velocity(d) = MutantCol.Velocity(d);
                    else
                        TrialCol.Velocity(d) = wCol.Velocity(d);
                    end
                end

                % Constraints are applied to the velocities            
                VelMax = +Alpha.*((up-down)./up); VelMin = -VelMax;                     
                TrialCol.Velocity = min(max(TrialCol.Velocity,VelMin),VelMax);

                TrialCol.Position = wCol.Position + TrialCol.Velocity;

                for flg=1:nVar
                    if TrialCol.Position(flg) < down | TrialCol.Position(flg) > up
                        TrialCol.Velocity(flg) = -TrialCol.Velocity(flg);              
                    end
                end
                for flg=1:nVar
                        TrialCol.Position(flg) = min(up,max(down,TrialCol.Position(flg))); % Bound the new location
                end                    

                TrialCol.Cost = CostFunction(TrialCol.Position);
                FEs = FEs + 1;
%                 % Selection
%                 if (TrialCol.Cost*(norm(emp(k).Imp.Position-TrialCol.Position))) < (wCol.Cost*(norm(emp(k).Imp.Position-wCol.Position)))
%                     emp(k).Col(indexwCol).Position = TrialCol.Position;                    
%                     emp(k).Col(indexwCol).Velocity = TrialCol.Velocity;
%                     emp(k).Col(indexwCol).Cost = TrialCol.Cost;
%                     emp(k).Col(indexwCol).Power = TrialCol.Power;
%                 end

                % Update personal best of the selected colony
%                 if emp(k).Col(indexwCol).Power >= emp(k).Col(indexwCol).Best.Power
%                     emp(k).Col(indexwCol).Best.Position = emp(k).Col(indexwCol).Position;
%                     emp(k).Col(indexwCol).Best.Velocity = emp(k).Col(indexwCol).Velocity;
%                     emp(k).Col(indexwCol).Best.Cost = emp(k).Col(indexwCol).Cost;
%                     emp(k).Col(indexwCol).Best.Power = emp(k).Col(indexwCol).Power;
%                 end

%                 % Update the Global Best imperialist using colonies
%                 if emp(k).Col(indexwCol).Power >= GlobalBest.Power
%                     temp = GlobalBest;
%                     GlobalBest.Position = emp(k).Col(indexwCol).Position;
%                     GlobalBest.Velocity = emp(k).Col(indexwCol).Velocity;
%                     GlobalBest.Cost = emp(k).Col(indexwCol).Cost;
%                     GlobalBest.Power = emp(k).Col(indexwCol).Power;
%                     emp(k).Col(indexwCol).Position = temp.Position;
%                     emp(k).Col(indexwCol).Velocity = temp.Velocity;
%                     emp(k).Col(indexwCol).Cost = temp.Cost;
%                     emp(k).Col(indexwCol).Power = temp.Power;
%                 end

%                 % Exchange a colony with its local imperialist if the colony has more power than its local imperialist
%                 if emp(k).Col(indexwCol).Power > emp(k).Imp.Power
%                     [emp(k).Imp, emp(k).Col(indexwCol)] = deal(emp(k).Col(indexwCol), emp(k).Imp);                
%                 end
        end
    end
end

%% Bring out all colonies and imperialists from their relevant empires and store in archive F ...
country=repmat(empty_country,1,1);
for k = 1:numel(emp)
    A = emp(k).Imp;
    country = [country; A]; %#ok<*AGROW>
    if numel(emp(k).Col) > 0
        for col = 1:numel(emp(k).Col)
            B = emp(k).Col(col);
            country = [country; B];
        end
    end
end
tempGL.Best.Position = tempGL.Position;
tempGL.Best.Velocity = tempGL.Velocity;
tempGL.Best.Cost = tempGL.Cost;
tempGL.Best.Power = tempGL.Power;
tempGL.Best.Rank = tempGL.Rank;
tempGL.Best.DominationSet = tempGL.DominationSet;
tempGL.Best.DominatedCount = tempGL.DominatedCount;
tempGL.Best.PS = tempGL.PS;
country = [country; tempGL];

country(numel(country)+1).Position = TrialCol.Position;
country(end).Velocity = TrialCol.Velocity;
country(end).Cost = TrialCol.Cost;
TrialCol.Power = 0; country(end).Power = TrialCol.Power;
TrialCol.Rank = 0; country(end).Rank = TrialCol.Rank;
TrialCol.DominationSet = 0; country(end).DominationSet = TrialCol.DominationSet;
TrialCol.DominatedCount = 0; country(end).DominatedCount = TrialCol.DominatedCount;
TrialCol.PS = 0; country(end).PS = TrialCol.PS;

TrialCol.Best.Position = TrialCol.Position;
TrialCol.Best.Velocity = TrialCol.Velocity;
TrialCol.Best.Cost = TrialCol.Cost;
TrialCol.Best.Power = TrialCol.Power;
TrialCol.Best.Rank = TrialCol.Rank;
TrialCol.Best.DominationSet = 0;
TrialCol.Best.DominatedCount = 0;
TrialCol.Best.PS = 0;
country(end).Best = TrialCol.Best;

country(1, :) = [];

%% Fast Non-dominated sorting, Penalized Moment Density estimating associated with Ranking, Fitness (Power) estimating, Sorting and Sub-population forming ...
[country, F] = NonDominatedSorting(country);
% Calculate Density Estimation in objective and decision spaces
country = PS(country,F);
% Calculate Fitness (Power)
country = PowerEstimator(country);
% Sort countries
[~,index] = sort([country.Power],'descend');
country=country(index);

%%  Form imperialists and colonies    

% Sort countries
[~,index] = sort([country.Power],'descend');
country=country(index);
% Assign of Colonies and Imperialists
imp=country(1:nEmp);
col=country(nEmp+1:nPop);
empty_empire.Imp=[];
empty_empire.Col=repmat(empty_country,0,1);
empty_empire.nCol=0;
emp=repmat(empty_empire,nEmp,1);

% Assign Imperialists
    for k=1:nEmp
        emp(k).Imp=imp(k);
    end
        
% Assign Colonies
Np = [imp.Power]./max([imp.Power]);
SNp = sum(Np);
Npp = Np./SNp;

for Num=1:nCol
     k = RouletteWheelSelection(Npp);
     check = isempty(emp(k).Col);
     if check == 1
         emp(k).Col = col(Num);
     else
        emp(k).Col = [emp(k).Col; col(Num)];
     end
     emp(k).nCol = emp(k).nCol + 1;                 
end

index = find([emp.nCol]==0);
if index>=1
    for tedademp=1:numel(index)
        [quan, ind]=max([emp.nCol]);
        randcol = randi(quan);
        emp(index(tedademp)).Col = emp(ind).Col(randcol);
        emp(index(tedademp)).nCol = 1;
        emp(ind).Col(randcol) = [];
        emp(ind).nCol = emp(ind).nCol - 1;             
    end  
end

%% Personal best, Global best, and Archive best Update, Exchange a colony with its imperialist if colony’s fitness (power) is more ...    
for k = 1:numel(emp)
    for col = 1:numel(emp(k).Col)
            % Update personal best of the colonies
            if emp(k).Col(col).Power >= emp(k).Col(col).Best.Power
                emp(k).Col(col).Best.Position = emp(k).Col(col).Position;
                emp(k).Col(col).Best.Velocity = emp(k).Col(col).Velocity;
                emp(k).Col(col).Best.Cost = emp(k).Col(col).Cost;                    
                emp(k).Col(col).Best.Power = emp(k).Col(col).Power;
                emp(k).Col(col).Best.Rank = emp(k).Col(col).Rank;
                emp(k).Col(col).Best.DominationSet = emp(k).Col(col).DominationSet;
                emp(k).Col(col).Best.DominatedCount = emp(k).Col(col).DominatedCount;                
                emp(k).Col(col).Best.PS = emp(k).Col(col).PS;
            end
            % Update the Global Best imperialist using colonies
            if emp(k).Col(col).Power >= GlobalBest.Power
                temp = GlobalBest;
                GlobalBest.Position = emp(k).Col(col).Position;
                GlobalBest.Velocity = emp(k).Col(col).Velocity;
                GlobalBest.Cost = emp(k).Col(col).Cost;                    
                GlobalBest.Power = emp(k).Col(col).Power;
                GlobalBest.Rank = emp(k).Col(col).Rank;            
                GlobalBest.DominationSet = emp(k).Col(col).DominationSet;  
                GlobalBest.DominatedCount = emp(k).Col(col).DominatedCount;                            
                GlobalBest.PS = emp(k).Col(col).PS;                                                                        

                emp(k).Col(col).Position = temp.Position;
                emp(k).Col(col).Velocity = temp.Velocity;
                emp(k).Col(col).Cost = temp.Cost;                    
                emp(k).Col(col).Power = temp.Power;
                emp(k).Col(col).Rank = temp.Rank;
                emp(k).Col(col).DominationSet = temp.DominationSet;
                emp(k).Col(col).DominatedCount = temp.DominatedCount;                
                emp(k).Col(col).PS = temp.PS;                                                                                                    
            end

            % Exchange a colony with its local imperialist if the colony has more power than its local imperialist
            if emp(k).Col(col).Power > emp(k).Imp.Power
                [emp(k).Imp, emp(k).Col(col)] = deal(emp(k).Col(col), emp(k).Imp);                
            end

            if emp(k).Imp.Power >= emp(k).Imp.Best.Power
                emp(k).Imp.Best.Position = emp(k).Imp.Position;
                emp(k).Imp.Best.Velocity = emp(k).Imp.Velocity;
                emp(k).Imp.Best.Cost = emp(k).Imp.Cost;            
                emp(k).Imp.Best.Power = emp(k).Imp.Power;
                emp(k).Imp.Best.Rank = emp(k).Imp.Rank;
                emp(k).Imp.Best.DominationSet = emp(k).Imp.DominationSet;
                emp(k).Imp.Best.DominatedCount = emp(k).Imp.DominatedCount;
                emp(k).Imp.Best.PS = emp(k).Imp.PS;                                    
            end

            % Update the Global Best imperialist using imperialist
            if emp(k).Imp.Power >= GlobalBest.Power
                temp = GlobalBest;
                GlobalBest.Position = emp(k).Imp.Position;
                GlobalBest.Velocity = emp(k).Imp.Velocity;
                GlobalBest.Cost = emp(k).Imp.Cost;            
                GlobalBest.Power = emp(k).Imp.Power;
                GlobalBest.Rank = emp(k).Imp.Rank;            
                GlobalBest.DominationSet = emp(k).Imp.DominationSet;            
                GlobalBest.DominatedCount = emp(k).Imp.DominatedCount;            
                GlobalBest.PS = emp(k).Imp.PS;            

                emp(k).Imp.Position = temp.Position;
                emp(k).Imp.Velocity = temp.Velocity;
                emp(k).Imp.Cost = temp.Cost;            
                emp(k).Imp.Power = temp.Power;
                emp(k).Imp.Rank = temp.Rank;
                emp(k).Imp.DominationSet = temp.DominationSet;
                emp(k).Imp.DominatedCount = temp.DominatedCount;                
                emp(k).Imp.PS = temp.PS;            
            end                                                
    end
end

%% Save the best 1st ranked solutions in the Archive
for k = 1:numel(emp)
    if numel(emp(k).Col)>0
        for col = 1 : numel(emp(k).Col)
             Archive = [Archive; emp(k).Col(col)];
        end % Colonies and imperialists move toward their best positions and
        % will be able to effectively update the Archive
    end
   Archive = [Archive; emp(k).Imp];
end

   Archive = [Archive; TrialCol];
   Archive = [Archive; tempGL];

 if it == 1
     Archive = Archive(2:end);
 end
%% Save Global Best results

%BestPower(it) = GlobalBest.Power;
BestPower(it) = max([Archive.Power]); 
    if FEs >= MaxFEs
        break;
    end
end
%% Save FEs in each Run
FEsRuns(nRun) = FEs;

% Fast Non-Dominated Sorting on Archive
[Archive, F] = NonDominatedSorting(Archive);

% Calculate Density Estimation in objective space
Archive = PS(Archive,F);

% Calculate Fitness (Power)
Archive = PowerEstimator(Archive);

Temp = [Archive.Rank];
Archive = Archive(find(Temp == 1)); %#ok<*FNDSB>

% Sort the 1st ranked solutions in the Archive based on their power values
[~, index] = sort([Archive.Power], 'descend');
Archive = Archive(index);

% % Set the size of Archive as nPop/2
% if numel(Archive) > nPop/2
%     Archive = Archive(1:nPop/2);
% end

%% Save IGD in each Run
PF = [Archive.Cost];
f1tf = down:1/size(PF,2):up;
truePF= 1-sqrt(f1tf);
if size(truePF,2) > size(PF,2)
    truePF = truePF(1:end,1: size(PF,2));
elseif size(PF,2) > size(truePF,2)
    PF = PF(1:end,1: size(truePF,2));
end
IGDRuns(nRun) = IGD(truePF, PF);
%% Save NHV in each Run
NHVRuns(nRun) = NHV(truePF, PF);

end

%% Calculate IGD results
BestIGDResults = min(IGDRuns);
MeanIGDResults = mean(IGDRuns);
MedianIGDResults = median(IGDRuns);
SDIGDResults = std(IGDRuns);

%% Show IGD results
disp('With this adjustment, the results are :');
disp([ ' Best IGD = '  num2str(BestIGDResults)]);
disp([ ' Mean IGD = '  num2str(MeanIGDResults)]);
disp([ ' Median IGD = '  num2str(MedianIGDResults)]);
disp([ ' SD IGD = '  num2str(SDIGDResults)]);

%% Save IGD results
save('BestIGDResults.mat','BestIGDResults');
save('MeanIGDResults.mat','MeanIGDResults');
save('MedianIGDResults.mat','MedianIGDResults');
save('SDIGDResults.mat','SDIGDResults');

%% Calculate NHV results
BestNHVResults = max(NHVRuns);
MeanNHVResults = mean(NHVRuns);
MedianNHVResults = median(NHVRuns);
SDNHVResults = std(NHVRuns);

%% Show NHV results
disp([ ' Best NHV = '  num2str(BestNHVResults)]);
disp([ ' Mean NHV = '  num2str(MeanNHVResults)]);
disp([ ' Median NHV = '  num2str(MedianNHVResults)]);
disp([ ' SD NHV = '  num2str(SDNHVResults)]);

%% Save NHV results
save('BestNHVResults.mat','BestNHVResults');
save('MeanNHVResults.mat','MeanNHVResults');
save('MedianNHVResults.mat','MedianNHVResults');
save('SDNHVResults.mat','SDNHVResults');

%% For convergence speed analysis, we calculated the total FEs considering the stopping criterion of "global best's power tolerance<=10^-2"
%BestFEsRunha = min(FEsRuns);
%MeanFEsRunha = mean(FEsRuns);