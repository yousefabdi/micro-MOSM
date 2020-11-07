
function Rep=CrowdingDistance(Rep)    
    
    % Sort Repository
    tmp=[Rep.Cost];
    [~,ind]=sort(tmp(1,:));
    Rep=Rep(ind);
    
    costs=[Rep.Cost]';
    [~,I,~] = unique(costs, 'rows', 'first'); 
    Rep=Rep(I);
    CrowdDis = zeros(1,numel(Rep));
    if numel(Rep)>=2
        costs=[Rep.Cost]';
        Fmax = max(costs,[],1);
        Fmin = min(costs,[],1);
        for i=1:size(costs,2)
           [~,rank] = sortrows(costs(:,i));
            CrowdDis(rank(1))   = inf;
            CrowdDis(rank(end)) = inf; 
            for j = 2 : size(costs,1)-1
                if (Fmax(i)-Fmin(i))>0
                    CrowdDis(rank(j)) = CrowdDis(rank(j))+(costs(rank(j+1),i)-costs(rank(j-1),i))/(Fmax(i)-Fmin(i));
                else
                    CrowdDis(rank(j)) = CrowdDis(rank(j))+(costs(rank(j+1),i)-costs(rank(j-1),i));
                end
            end
        end
        for j = 1 : size(costs,1)
            Rep(j).CD=CrowdDis(j);
        end
    elseif numel(Rep)==1
        Rep(1).CD=inf;
    end
 
end