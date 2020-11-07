
function Rep=MakeRepository(Structure, Rep)
    
    global SM_Settings;
    nRep = SM_Settings.nRep;
  
    Temp_Rep=[];
                                       
    Temp_Rep=[Temp_Rep
    Structure(~[Structure.IsDominated])];                

    Temp_Rep=[Temp_Rep
        Rep];
    Costs=[Temp_Rep.Cost]; 
    for i=1:numel(Temp_Rep)
        if numel(Temp_Rep)<2
            break;
        end
        y=Temp_Rep(i).Cost;
        Costs=[Temp_Rep.Cost];
        res=find(all(Costs<=y) & any(Costs<y));
        if ~isempty(res)
                Temp_Rep(i).IsDominated=true;
        end        
        
    end
    Rep=Temp_Rep(~[Temp_Rep.IsDominated]);
    
    Rep=CrowdingDistance(Rep);
    CDs=[Rep.CD];
    if size(Rep,1)>nRep  
        [~,rank] = sort(CDs,'descend');
        Rep = Rep(rank(1:min(nRep,end)));
        
        Rep=CrowdingDistance(Rep);        
    end            

end