
function indices=RouletteWheelSelection(P, n)

    r=rand;
    
    C=cumsum(P);
    indices=[];
    for j=1:n
        r=rand;
        i=find(r<=C,1,'first');              
        indices = [indices i];
        %C(i)=-1; 
    end
        
end