
function PlotCosts2(rep)
    
    p_costs=[];
    nGroup=4;
%     for i=1:nGroup
%         for j=1:numel(Structure(i).Group)
%            p_costs=[p_costs  Structure(i).Group(j).Cost];
%         end
%     end
    figure(1);
    
    rep_costs=rep;

    plot(rep_costs(1,:),rep_costs(2,:),'bo');
    hold on;
    %plot(p_costs(1,:),p_costs(2,:),'r+');
    %hold on;
    xlabel('F_1');
    ylabel('F_2');
    
    hold off;

end