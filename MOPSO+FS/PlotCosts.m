%
% Copyright (c) 2019, Qianying Liu
% All rights reserved. Please read the "license.txt" for license terms.
% Project Title: 
% Contact Info: 
%已改第一版，rep输出按特征个数排序
%

function [rep_costs]=PlotCosts(rep,nRep)  
    rep_costs=[rep.Cost];
    [featureNum,index]=sort(rep_costs(1,:),2);
    error=rep_costs(2,index);
    rep_costs=[featureNum;error]
    plot(featureNum,error,'r*');
    
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    title('MOPSO+FS');
    
    grid on;  
    hold off;
    
    rep_costs=zeros(2,nRep);
    rep_costs(:,1:numel(rep))=[featureNum;error];
end