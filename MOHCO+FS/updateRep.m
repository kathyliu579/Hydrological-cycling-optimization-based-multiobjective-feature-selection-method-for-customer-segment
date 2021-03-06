function [rep,pop] = updateRep(rep,pop,nGrid,alpha,nRep)
     % Determine Domination
    pop=DetermineDomination(pop);
    rep=[rep;pop(~[pop.IsDominated])];  % Add Non-Dominated Particles to REPOSITORY
    rep=DetermineDomination(rep);  % Determine Domination of New Resository Members
    rep=rep(~[rep.IsDominated]);  % Keep only Non-Dminated Memebrs in the Repository
    %delete the repetitive indvidual
    n=numel(rep);index=[];
    for i=1:n-1
        for j=i+1:n
            if (rep(i).Cost==rep(j).Cost)
                index=[index,j];
            end
        end
    end
    rep(index)=[];
    
    % Update Grid
    Grid=CreateGrid(rep,nGrid,alpha); 
    
    % Update Grid Indices
    for j=1:numel(rep)
        rep(j)=FindGridIndex(rep(j),Grid);
    end
    
    % Check if Repository is Full
    if numel(rep)>nRep
        Extra=numel(rep)-nRep;
        for e=1:Extra
            rep=DeleteOneRepMemebr(rep,gamma);
        end
    end
end

