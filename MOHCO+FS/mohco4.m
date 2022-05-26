% Copyright (c) 2019, Qianying Liu
% All rights reserved. Please read the "license.txt" for license terms.
% Project Title:
% Contact Info:
%已改第一版，贪婪规则，向rep学习，基于distance的蒸发已加
%流动 only learn from rep（加上pbest 不好）
%boundary process(0.5p 扰动)
%% repetited
clc;
clear;
close all;
%% Problem Definition
%Ionosphere_alter,WBCD_alter,Wine_alter,Zoo_alter
name={'Ionosphere_alter','WBCD_alter','Wine_alter','Zoo_alter','Australian','Default','German'};
for pronum=6
    data=LoadData(name(pronum));
    CostFunction=@(s) FeatureSelectionCost(s,data);      % Cost Function
    nVar=data.nx;             % Number of Decision Variables
    VarSize=[1 nVar];   % Size of Decision Variables Matrix
    VarMin=0;          % Lower Bound of Variables
    VarMax=1;          % Upper Bound of Variables
    %% MOHCO Parameters
    maxEvaluation=20000;           % Maximum Number of Iterations
    run=25;
    nPop=50;            % Population Size
    nRep=data.nx;            % Repository Size
    maxFlowTimes=3;     % maxFlowTimes
    par2=0.1;           % evepotature pobility
    minValue=0;
    hco4=zeros(2*run,nRep);
    c1=1;
    c2=0;
    %% grid Parameters
    nGrid=7;            % Number of Grids per Dimension
    alpha=0.1;          % Inflation Rate
    beta=2;             % Leader Selection Pressure
    gamma=2;            % Deletion Selection Pressure
    
    for runNow=1:run
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
        while(terminate(evaluationNow,maxEvaluation))   %终止函数
            %% 流动算子
            for i=1:nPop
                nrep=numel(rep);
                % 计算流动后的位置及适应度值
                NewSol.Position=pop(i).Position+c1*rand(1,nVar).*(rep(randi(nrep)).Position-pop(i).Position)+c2*rand(1,nVar).*(pop(i).Best.Position-pop(i).Position);   %向所学习的流动方向的粒子进行学习
                NewSol.Cost=CostFunction(NewSol.Position) ;           %重新计算流动后的适应值
                evaluationNow=evaluationNow+1;                       %评价后的评估次数+1
                %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
                indexmin=NewSol.Position<VarMin;
                indexmax=NewSol.Position>VarMax;
                
                if rand<0.5
                    NewSol.Position(indexmin) = VarMin;
                    NewSol.Position(indexmax) = VarMax;
                else
                    mi=sum(indexmin); ma=sum(indexmax);
                    a=VarMin+rand(1,mi)*0.2;   b=VarMax-rand(1,ma)*0.2;
                    NewSol.Position(indexmin) =a;
                    NewSol.Position(indexmax)= b;
                end
                if NewSol.Position ==0
                    index=randi(nVar);
                    NewSol.Position(index)=rand;
                end
                flowtimes=0;
                while flowtimes<maxFlowTimes
                    if NewSol.Cost<pop(i).Cost                       %如果新生成的适应值小于粒子原本的适应值
                        pop(i).Position=NewSol.Position;             %则该粒子的位置被替换
                        pop(i).Cost=NewSol.Cost;                     %该例子的适应度值被替换
                        NewSol.Position=pop(i).Position+c1*rand(1,nVar).*(rep(randi(nrep)).Position-pop(i).Position)+c2*rand(1,nVar).*(pop(i).Best.Position-pop(i).Position);  %向所学习的流动方向的粒子进行学习
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
                            
                            NewSol.Position=pop(i).Position+rand(1,nVar).*(rep(randi(nrep)).Position-pop(i).Position);   %向所学习的流动方向的粒子进行学习
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
                %    boundary process
                indexmin=NewSol.Position<VarMin;
                indexmax=NewSol.Position>VarMax;
                
                if rand<0.5
                    NewSol.Position(indexmin) = VarMin;
                    NewSol.Position(indexmax) = VarMax;
                else
                    mi=sum(indexmin); ma=sum(indexmax);
                    a=VarMin+rand(1,mi)*0.2;   b=VarMax-rand(1,ma)*0.2;
                    NewSol.Position(indexmin) =a;
                    NewSol.Position(indexmax)= b;
                end
                if NewSol.Position ==0
                    index=randi(nVar);
                    NewSol.Position(index)=rand;
                end
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
            pEva=createPossibility(pop,nPop);
            for i=1:nPop
                if rand()<pEva(i)
                    if rand()<0.1
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
                        indexmin=newIndividual<VarMin;
                        indexmax=newIndividual>VarMax;
                        
                        if rand<0.5
                            newIndividual(indexmin) = VarMin;
                            newIndividual(indexmax) = VarMax;
                        else
                            mi=sum(indexmin); ma=sum(indexmax);
                            a=VarMin+rand(1,mi)*0.2;   b=VarMax-rand(1,ma)*0.2;
                            newIndividual(indexmin) =a;
                            newIndividual(indexmax)= b;
                        end
                        if newIndividual ==0
                            index=randi(nVar);
                            newIndividual(index)=rand;
                        end
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
            [hco_costs]=PlotCosts(rep,nRep);
            % PlotCosts(pop,rep);
            pause(0.01);
            % Show Iteration Information
            disp(['Evaluation ' num2str(evaluationNow) ': Number of Rep Members = ' num2str(numel(rep))]);
            
        end
        [hco_costs]=PlotCosts(rep,nRep);
        hco4((runNow-1)*2+1:(runNow-1)*2+2,:)=hco_costs;
        if pronum==1
            save('hco4_Ionosp','hco4')
        elseif pronum==2
            save('hco4_WBCD','hco4')
        elseif pronum==3
            save('hco4_Wine','hco4')
        elseif pronum==4
            save('hco4_Zoo','hco4')
        elseif pronum==5
            save('hco4_Australian','hco4')
        elseif pronum==6
            save('hco4_Default','hco4')
        else
            save('hco4_German','hco4')
        end
    end
end

