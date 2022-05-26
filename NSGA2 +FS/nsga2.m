%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA120
% Project Title: Non-dominated Sorting Genetic Algorithm II (NSGA-II)
% Publisher: Yarpiz (www.yarpiz.com)
%
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
%
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%
clc;
clear;
close all;
%% Problem Definition
name={'Ionosphere_alter','WBCD_alter','Wine_alter','Zoo_alter','Australian','Default','German'};
for pronum=7
    data=LoadData(name(pronum));
    
    CostFunction=@(s) FeatureSelectionCost(s,data);      % Cost Function
    
    nVar=data.nx;             % Number of Decision Variables 在特征选择问题中就是数据集本身特征的总数
    
    VarSize=[1 nVar];   % Size of Decision Variables Matrix
    
    VarMin=0;          % Lower Bound of Variables
    VarMax= 1;          % Upper Bound of Variables
    
    % Number of Objective Functions
    nObj=numel(CostFunction(randi([0 1],VarSize)));
    %% NSGA-II Parameters
    run=4;
    MaxIt=100;      % Maximum Number of Iterations
    nsga=[];
    
    nRep=data.nx;
    nPop=50;        % Population Size
    pCrossover=0.5;                         % Crossover Percentage
    nCrossover=2*round(pCrossover*nPop/2);  % Number of Parnets (Offsprings)
    
    pMutation=0.5;                          % Mutation Percentage
    nMutation=round(pMutation*nPop);        % Number of Mutants
    mu=0.02;                    % Mutation Rate
    sigma=0.1*(VarMax-VarMin);  % Mutation Step Size
    
    for runNow=run
        %% Initialization
        
        empty_individual.Position=[];
        empty_individual.Cost=[];
        empty_individual.Rank=[];
        empty_individual.DominationSet=[];
        empty_individual.DominatedCount=[];
        empty_individual.CrowdingDistance=[];
        
        pop=repmat(empty_individual,nPop,1);
        
        for i=1:nPop
            
            if i~=1
                pop(i).Position=randi([0 1],VarSize);
            else
                pop(i).Position=ones(VarSize);
            end
            
            [pop(i).Cost, pop(i).Out]=CostFunction(pop(i).Position);
            
        end
        
        % Non-Dominated Sorting
        [pop, F]=NonDominatedSorting(pop);
        
        % Calculate Crowding Distance
        pop=CalcCrowdingDistance(pop,F);
        
        % Sort Population
        [pop, F]=SortPopulation(pop);
        
        %% NSGA-II Main Loop
        for it=1:MaxIt
            % Crossover
            popc=repmat(empty_individual,nCrossover/2,2);
            for k=1:nCrossover/2
                
                i1=randi([1 nPop]);
                p1=pop(i1);
                
                i2=randi([1 nPop]);
                p2=pop(i2);
                
                [popc(k,1).Position, popc(k,2).Position]=Crossover(p1.Position,p2.Position);
                
                [popc(k,1).Cost, popc(k,1).Out]=CostFunction(popc(k,1).Position);
                [popc(k,2).Cost, popc(k,2).Out]=CostFunction(popc(k,2).Position);
                
                
            end
            popc=popc(:);
            
            % Mutation
            popm=repmat(empty_individual,nMutation,1);
            for k=1:nMutation
                
                i=randi([1 nPop]);
                p=pop(i);
                
                popm(k).Position=Mutate(p.Position,mu,sigma);
                [popm(k).Cost, popm(k).Out]=CostFunction(popm(k).Position);
            end
            
            % Merge
            pop=[pop
                popc
                popm]; %#ok
            
            % Non-Dominated Sorting
            [pop, F]=NonDominatedSorting(pop);
            
            % Calculate Crowding Distance
            pop=CalcCrowdingDistance(pop,F);
            
            % Sort Population
            [pop, F]=SortPopulation(pop); %#ok
            
            % Truncate
            pop=pop(1:nPop);
            
            % Non-Dominated Sorting
            [pop, F]=NonDominatedSorting(pop);
            
            % Calculate Crowding Distance
            pop=CalcCrowdingDistance(pop,F);
            
            % Sort Population
            [pop, F]=SortPopulation(pop);
            
            % Store F1
            F1=pop(F{1});
            F1=GetUniqueMembers(F1);
            % Show Iteration Information
            disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);
            
            % Plot F1 Costs
            figure(1);
            [ga_costs]=PlotCosts(F1,nRep);
            
            pause(0.01);
        end
        nsga((runNow-1)*2+1:(runNow-1)*2+2,:)=ga_costs;
        
        if pronum==1
            save('nsga_Ionosp_1run','nsga')
        elseif pronum==2
            save('nsga_WBCD_1run','nsga')
        elseif pronum==3
            save('nsga_Wine','nsga')
        elseif pronum==4
            save('nsga_Zoo','nsga')
        elseif pronum==5
            save('nsga_Australian','nsga')
        elseif pronum==6
            save('nsga_Default','nsga')
        else
            save('nsga_German_new4','nsga')
        end
    end
end


%% Results

