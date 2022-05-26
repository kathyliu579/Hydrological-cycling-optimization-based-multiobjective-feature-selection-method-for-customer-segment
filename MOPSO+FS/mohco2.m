%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA121
% Project Title: Multi-Objective Particle Swarm Optimization (MOPSO)
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

data=LoadData();

CostFunction=@(s) FeatureSelectionCost(s,data);      % Cost Function

nVar=data.nx;             % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=0;          % Lower Bound of Variables
VarMax=1;          % Upper Bound of Variables


%% MOPSO Parameters

maxEvaluation=10000;           % Maximum Number of Iterations

nPop=50;            % Population Size

nRep=15;            % Repository Size

maxFlowTimes=3;     % maxFlowTimes
par2=0.1;           % evepotature pobility
minValue=0;

nGrid=7;            % Number of Grids per Dimension
alpha=0.1;          % Inflation Rate

beta=2;             % Leader Selection Pressure
gamma=2;            % Deletion Selection Pressure

% mu=0.1;             % Mutation Rate

%% Initialization

empty_particle.Position=[];
empty_particle.Velocity=[];
empty_particle.Cost=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
empty_particle.IsDominated=[];
empty_particle.GridIndex=[];
empty_particle.GridSubIndex=[];

pop=repmat(empty_particle,nPop,1);

% Determine ImproveTrial/evaluation/evaluationNow
improveTrial=zeros(1,nPop);      improveTrial2=improveTrial;                        %提高试验
evaluation=0;                                            %%评估次数归零
evaluationNow=0;                                         %%记录评价后的评估次数

for i=1:nPop
    
    if i~=1
        pop(i).Position=randi([0 1],VarSize);
    else
        pop(i).Position=ones(VarSize);
    end
    
    pop(i).Velocity=zeros(VarSize);
    
    pop(i).Cost=CostFunction(pop(i).Position);
    
    % Update Personal Best
    pop(i).Best.Position=pop(i).Position;
    pop(i).Best.Cost=pop(i).Cost;
end

evaluationNow=evaluationNow+nPop;

% Determine Domination
pop=DetermineDomination(pop);

rep=pop(~[pop.IsDominated]);

Grid=CreateGrid(rep,nGrid,alpha);

for i=1:numel(rep)
    rep(i)=FindGridIndex(rep(i),Grid);
end


%% MOPSO Main Loop

%迭代次数=0
while(terminate(evaluationNow,maxEvaluation))   %终止函数
    
    %% 流动算子
    for i=1:nPop
        
        nrep=numel(rep);
        
        if pop(i).Cost==rep(randi(nrep)).Cost
            flowDirection=fix(rand*(nPop))+1;   %fix表示向下取整
            while flowDirection==i
                flowDirection=fix(rand*(nPop))+1;
            end
        else
            for j=1:nPop
                f1(j)=pop(j).Cost(1);
                f2(j)=pop(j).Cost(2);
            end
            betterInd=find(f1<=pop(j).Cost(1)&f2<=pop(j).Cost(2));
            
            if isempty(betterInd)==0 %exist better one
                flowDirection=betterInd(fix(rand*(length(betterInd)))+1);
            else
                flowDirection=fix(rand*(nPop))+1;
            end
        end
        % 计算流动后的位置及适应度值
        
        NewSol.Position=pop(i).Position+0.5*rand(1,nVar).*(pop(flowDirection).Position-pop(i).Position)+0.5*rand(1,nVar).*(rep(randi(nrep)).Position-pop(i).Position);   %向所学习的流动方向的粒子进行学习
        NewSol.Cost=CostFunction(NewSol.Position) ;           %重新计算流动后的适应值
        evaluationNow=evaluationNow+1;                       %评价后的评估次数+1
        %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
        NewSol.Position = max(NewSol.Position, VarMin);
        NewSol.Position = min(NewSol.Position, VarMax);
        
        flowtimes=0;
        while flowtimes<maxFlowTimes
            if NewSol.Cost<pop(i).Cost                       %如果新生成的适应值小于粒子原本的适应值
                pop(i).Position=NewSol.Position;             %则该粒子的位置被替换
                pop(i).Cost=NewSol.Cost;                     %该例子的适应度值被替换
                
                NewSol.Position=pop(i).Position+0.5*rand(1,nVar).*(pop(flowDirection).Position-pop(i).Position)+rand(1,nVar).*(rep(randi(nrep)).Position-pop(i).Position);   %向所学习的流动方向的粒子进行学习
                NewSol.Cost=CostFunction(pop(i).Position);            %重新计算流动后的适应值
                evaluationNow=evaluationNow+1;
                flowtimes=flowtimes+1;
            elseif NewSol.Cost>pop(i).Cost
                improveTrial(i)=improveTrial(i)+1;          %改进试验+1
                flowtimes=maxFlowTimes;
            else
                if rand<0.5
                    pop(i).Position=NewSol.Position;             %则该粒子的位置被替换
                    pop(i).Cost=NewSol.Cost;                     %该例子的适应度值被替换
                    
                    NewSol.Position=pop(i).Position+0.5*rand(1,nVar).*(pop(flowDirection).Position-pop(i).Position)+rand(1,nVar).*(rep(randi(nrep)).Position-pop(i).Position);   %向所学习的流动方向的粒子进行学习
                    NewSol.Cost=CostFunction(pop(i).Position);            %重新计算流动后的适应值
                    evaluationNow=evaluationNow+1;
                    flowtimes=flowtimes+1;
                else
                    improveTrial(i)=improveTrial(i)+1;
                end
            end
        end
        
    end
    %%     流动后更新最优位置及最优值
    %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
    pop(i).Position = max(pop(i).Position, VarMin);
    pop(i).Position = min(pop(i).Position, VarMax);
    
    if Dominates(pop(i),pop(i).Best)
        pop(i).Best.Position=pop(i).Position;
        pop(i).Best.Cost=pop(i).Cost;
        
    elseif Dominates(pop(i).Best,pop(i))
        % Do Nothing
        
    else
        if rand<0.5
            pop(i).Best.Position=pop(i).Position;
            pop(i).Best.Cost=pop(i).Cost;
        end
    end
    %%
    %update rep
    [rep,pop]=updateRep(rep,pop,nGrid,alpha,nRep);
    
    if(evaluationNow>maxEvaluation)
        break;
        %本次运行结束
    end
    
    %%     渗透操作
    for i=1:nPop
        %/*随机选择的解A用于产生对于粒子i的突变解
        %/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        neighbour=fix(rand*(nPop))+1;
        
        %/*Randomly selected solution must be different from the solution i*/
        while(neighbour==i)
            neighbour=fix(rand*(nPop))+1;
        end
        %%%%%每次所有的维度都变化%%%%%%
        
        NewSol.Position=pop(i).Position+(pop(i).Position-pop(neighbour).Position).*(rand(1,nVar)-0.5)*0.5;
        %             newIndividual=PositionPop(i,:)+(PositionPop(i,:)-PositionPop(neighbour,:)).*((rand(1,Dim)-0.5)*2);  %在选定的领域方向，在【-1，1】之间的噪音下学习渗透
        
        %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
        NewSol.Position = max(NewSol.Position, VarMin);
        NewSol.Position = min(NewSol.Position, VarMax);
        
        %evaluate new solution
        NewSol.Cost=CostFunction(NewSol.Position);
        evaluationNow=evaluationNow+1;                    %%适应度评估次数+1
        if(evaluationNow>maxEvaluation)
            break;
            %本次运行结束
        end
        
        % /*a greedy selection is applied between the current solution i and its mutant*/
        if (NewSol.Cost<pop(i).Cost)
            pop(i).Position=NewSol.Position;
            pop(i).Cost=NewSol.Cost;
        elseif (NewSol.Cost>pop(i).Cost)
            improveTrial2(i)=improveTrial2(i)+1;
        else
            if rand<0.5
                pop(i).Position=NewSol.Position;             %则该粒子的位置被替换
                pop(i).Cost=NewSol.Cost;                     %该例子的适应度值被替换
            else
                improveTrial2(i)=improveTrial2(i)+1;
            end
        end
    end
    %% 渗透后更新最优位置及最优值
    if Dominates(pop(i),pop(i).Best)
        pop(i).Best.Position=pop(i).Position;
        pop(i).Best.Cost=pop(i).Cost;
        
    elseif Dominates(pop(i).Best,pop(i))
        % Do Nothing
        
    else
        if rand<0.5
            pop(i).Best.Position=pop(i).Position;
            pop(i).Best.Cost=pop(i).Cost;
        end
    end
    %%% end of 渗透
    %%
    %update rep
    [rep,pop]=updateRep(rep,pop,nGrid,alpha,nRep);
    
    %% 蒸发和降水
    %随机蒸发？地势高蒸发？优劣各选部分蒸发？
    nrep=numel(rep);
    for i=1:nPop
        if rand()<par2
            if rand()<0.01
                pop(i).Position=VarMin+(VarMax-VarMin)*rand(1,nVar);  %%以0.2的概率随机蒸发和降雨
                pop(i).Cost=CostFunction(pop(i).Position);
                evaluationNow=evaluationNow+1;
                if(evaluationNow>maxEvaluation)
                    break;
                    %本次运行结束
                end
            else
                
                newIndividual=rep(randi(nrep)).Position;
                
                randDim=fix(rand*nVar)+1;
                gaussianPra=normrnd(0,1,[1,randDim])+1;   %产生正态分布随机数的函数normrnd-----用来产生高斯随机矩阵
                randOrder=randperm(nVar);
                for j=1:randDim
                    newIndividual(randOrder(j))=newIndividual(randOrder(j)).*gaussianPra(j);
                end
                
                %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
                newIndividual = max(newIndividual, VarMin);
                newIndividual = min(newIndividual, VarMax);
                
                pop(i).Position=newIndividual;  %%以0.2的概率随机蒸发和降雨
                pop(i).Cost=CostFunction(pop(i).Position);
            end
        end
        
    end
    %% 蒸发与降水后最优值与位置更新
    
    if Dominates(pop(i),pop(i).Best)
        pop(i).Best.Position=pop(i).Position;
        pop(i).Best.Cost=pop(i).Cost;
        
    elseif Dominates(pop(i).Best,pop(i))
        % Do Nothing
        
    else
        if rand<0.5
            pop(i).Best.Position=pop(i).Position;
            pop(i).Best.Cost=pop(i).Cost;
        end
    end
    
    %%
    %update rep
    [rep,pop]=updateRep(rep,pop,nGrid,alpha,nRep);
    
    %%
    % Plot Costs
    figure(1);
    PlotCosts(rep);
    %    PlotCosts(pop,rep);
    pause(0.01);
    
    % Show Iteration Information
    disp(['Evaluation ' num2str(evaluationNow) ': Number of Rep Members = ' num2str(numel(rep))]);
end
%%%%putout%%%%

% Plot Costs
figure(1);
PlotCosts(rep);
%    PlotCosts(pop,rep);
pause(0.01);

% Show Iteration Information
disp(['Evaluation ' num2str(evaluationNow) ': Number of Rep Members = ' num2str(numel(rep))]);
