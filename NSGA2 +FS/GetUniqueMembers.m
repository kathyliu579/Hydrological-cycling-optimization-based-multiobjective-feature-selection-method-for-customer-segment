%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPML122
% Project Title: Multi-Objective Feature selection using NSGA-II
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function UA=GetUniqueMembers(A)

    u=true(size(A));
    
    for i=2:numel(A)
        for j=1:i-1
            if IsSame(A(i),A(j))
                u(i)=false;
                break;
            end
        end
    end
    
    UA=A(u);

end

function e=IsSame(a,b)

    e=all(a.Position==b.Position);

end