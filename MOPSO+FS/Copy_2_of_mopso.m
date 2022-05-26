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
name={'Ionosphere_alter','WBCD_alter','Wine_alter','Zoo_alter','Australian','Default','German'};
for pronum=7
    data=LoadData(name(pronum));
    
CostFunction=@(s) FeatureSelectionCost(s,data);      % Cost Function

nVar=data.nx;             % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=0;          % Lower Bound of Variables
VarMax=1;          % Upper Bound of Variables
%% MOPSO Parameters
run=2;
MaxIt=100;           % Maximum Number of Iterations
pso=[];
nPop=50;            % Population Size

nRep=data.nx;            % Repository Size

w=0.5;              % Inertia Weight
wdamp=0.99;         % Intertia Weight Damping Rate
c1=1;               % Personal Learning Coefficient
c2=2;               % Global Learning Coefficient

nGrid=7;            % Number of Grids per Dimension
alpha=0.1;          % Inflation Rate

beta=2;             % Leader Selection Pressure
gamma=2;            % Deletion Selection Pressure

mu=0.1;             % Mutation Rate
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
    
    
    % Determine Domination
    pop=DetermineDomination(pop);
    
    rep=pop(~[pop.IsDominated]);
    
    Grid=CreateGrid(rep,nGrid,alpha);
    
    for i=1:numel(rep)
        rep(i)=FindGridIndex(rep(i),Grid);
    end
    
    
    %% MOPSO Main Loop
    
    for it=1:MaxIt
        
        for i=1:nPop
            
            leader=SelectLeader(rep,beta);
            
            pop(i).Velocity = w*pop(i).Velocity ...
                +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
                +c2*rand(VarSize).*(leader.Position-pop(i).Position);
            
            pop(i).Position = pop(i).Position + pop(i).Velocity;
            
            pop(i).Position = max(pop(i).Position, VarMin);
            pop(i).Position = min(pop(i).Position, VarMax);
            
            pop(i).Cost = CostFunction(pop(i).Position);
            
            % Apply Mutation
            pm=(1-(it-1)/(MaxIt-1))^(1/mu);
            if rand<pm
                NewSol.Position=Mutate(pop(i).Position,pm,VarMin,VarMax);
                NewSol.Cost=CostFunction(NewSol.Position);
                if Dominates(NewSol,pop(i))
                    pop(i).Position=NewSol.Position;
                    pop(i).Cost=NewSol.Cost;
                    
                elseif Dominates(pop(i),NewSol)
                    % Do Nothing
                    
                else
                    if rand<0.5
                        pop(i).Position=NewSol.Position;
                        pop(i).Cost=NewSol.Cost;
                    end
                end
            end
            
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
            
        end
                % Damping Inertia Weight
        w=w*wdamp;
        
        
        % Add Non-Dominated Particles to REPOSITORY
        rep=[rep
            pop(~[pop.IsDominated])]; %#ok
        %delete the repetitive indvidual
        n=numel(rep);index=[];
        for ii=1:n-1
            for jj=ii+1:n
                if (rep(ii).Cost==rep(jj).Cost)
                    index=[index,jj];
                end
            end
        end
        rep(index)=[];
        % Determine Domination of New Resository Members
        rep=DetermineDomination(rep);
        
        % Keep only Non-Dminated Memebrs in the Repository
        rep=rep(~[rep.IsDominated]);
        
        % Update Grid
        Grid=CreateGrid(rep,nGrid,alpha);
        
        % Update Grid Indices
        for i=1:numel(rep)
            rep(i)=FindGridIndex(rep(i),Grid);
        end
        
        % Check if Repository is Full
        if numel(rep)>nRep
            
            Extra=numel(rep)-nRep;
            for e=1:Extra
                rep=DeleteOneRepMemebr(rep,gamma);
            end
            
        end
        
        % Plot Costs
        figure(1);
        [pso_costs]=PlotCosts(rep,nRep);
        %    PlotCosts(pop,rep);
        pause(0.01);
        
        % Show Iteration Information
        disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);
   
    end
    pso((runNow-1)*2+1:(runNow-1)*2+2,:)=pso_costs;
    save('pso_Ionosp','pso')
       
    if pronum==1
        save('pso_Ionosp','pso')
    elseif pronum==2
        save('pso_WBCD','pso')
    elseif pronum==3
        save('pso_Wine','pso')
    elseif pronum==4
        save('pso_Zoo','pso')
    elseif pronum==5
        save('pso_Australian','pso')
    elseif pronum==6
        save('pso_Default_run1','pso')
    else
        save('pso_German_run2','pso')
    end
end
end
%% Resluts

