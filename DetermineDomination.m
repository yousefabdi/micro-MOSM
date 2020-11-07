
function Struct=DetermineDomination(Structure,sol_index)
    if sol_index==0
        nSol=numel(Structure);    

        for i=1:nSol        
            Structure(i).IsDominated=false;
        end

        for i=1:nSol

            Costs=[Structure.Cost];

            y=Structure(i).Cost;
            res=find(all(Costs<=y) & any(Costs<y));            
            if ~isempty(res)
                Structure(i).IsDominated=true;
            end
            
        end
        Struct=Structure;
    else
        for j=1:numel(Structure)
            Structure(j).IsDominated=false;
        end
        Costs=[Structure.Cost];

        y=Structure(sol_index).Cost;
        res=find(all(Costs<=y) & any(Costs<y));        
        if ~isempty(res)
            Structure(sol_index).IsDominated=true;
        end 
        Struct=Structure;
        
    end    

end