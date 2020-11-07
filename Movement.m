
function Structure=Movement(Structure, Rep, it)
    
    global ProblemSettings 
    CostFunction=ProblemSettings.CostFunction;
    VarSize=ProblemSettings.VarSize;
    VarMin=ProblemSettings.VarMin;
    MaxIt=ProblemSettings.MaxIt;
    nVar=ProblemSettings.nVar;
    VarMax=ProblemSettings.VarMax;  
    
    nSol=numel(Structure); 
    
    for i=1:nSol
        
        % Get the Group Style Vector
        Style=Structure(i).Style;

        % Compute Probability Vector for i'th Style Vector of Solution
        P=exp(-Style/sum(Style));
        P=1-P;
        P=P/sum(P);
        
        % Selection method parameters
        Wcd=0.5;
        Wsc=0.5;
        
        % Select one of the Groups based on the Probability Vector
        % and Roulette Wheel Selection Method       
        Selected_Style=RouletteWheelSelection(P,1);%
        
        %%%
            % The current Structure of Switch Cases Represents micro-MOSM-C1 
            % For micro-MOSM-C2: Case 1 -> Case 10
            %                    Case 10 -> Case 1
            %                    Case 2 -> Case 20
            %                    Case 20 -> Case 2
            
            % Note that the nStyle variable should be set to 3 in the SearchManager script
            % for micro-MOSM-C2
        %%%
        switch Selected_Style    % Which optimization methods should be run
                                 % at first they shoul be oprdered from 1                                                                 
                                    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 1 % Simulated binary crossover                
                costs=[Rep.Cost];
                cds=[Rep.CD];
                if all(cds==inf)
                    cds = ones(size(cds));
                else
                    cds(cds==inf) = 2*max(cds(cds~=inf));
                end
                s_costs=sum(costs,1);
                s_cds=sum(cds,1);
                P3=exp(-s_cds/sum(s_cds));
                P1=exp(-s_costs/sum(s_costs));
                P1=1-P1;
                P=Wcd*P3+Wsc*P1;
                P=P/sum(P);
                Selected_Sol=RouletteWheelSelection(P,1);%   

                Parent1=Structure(i);
                Parent2=Rep(Selected_Sol(1));

                beta = zeros(1,nVar);
                mu   = rand(1,nVar);
                [proC,disC,proM,disM] = deal(1,20,1,20);
                beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
                beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
                beta = beta.*(-1).^randi([0,1],1,nVar);
                beta(rand(1,nVar)<0.5) = 1;
                beta(repmat(rand(1,1)>proC,1,nVar)) = 1;
                Offspring = [(Parent1.Decision+Parent2.Decision)/2+beta.*(Parent2.Decision-Parent1.Decision)/2
                             (Parent2.Decision+Parent1.Decision)/2-beta.*(Parent2.Decision-Parent1.Decision)/2];
                      
                % Polynomial mutation
                Site  = rand(2*1,nVar) < proM/nVar;
                mu    = rand(2*1,nVar);
                temp  = Site & mu<=0.5;
                Vmin=VarMin;Vmax=VarMax;
                Vmin(2,:)=VarMin;Vmax(2,:)=VarMax;
                Offspring       = min(max(Offspring,VarMin),VarMax);
                Offspring(temp) = Offspring(temp)+(Vmax(temp)-Vmin(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                                  (1-(Offspring(temp)-Vmin(temp))./(Vmax(temp)-Vmin(temp))).^(disM+1)).^(1/(disM+1))-1);
                temp = Site & mu>0.5; 
                Offspring(temp) = Offspring(temp)+(Vmax(temp)-Vmin(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                                  (1-(Vmax(temp)-Offspring(temp))./(Vmax(temp)-Vmin(temp))).^(disM+1)).^(1/(disM+1)));
                                  
                Offspring = min(max(Offspring,VarMin),VarMax);   
                if rand<=0.5
                    Structure(i).Decision=Offspring(1,:);
                    Structure(i).PrevCost=Structure(i).Cost;
                    Structure(i).Cost=CostFunction(Offspring(1,:));
                else
                    Structure(i).Decision=Offspring(2,:);
                    Structure(i).PrevCost=Structure(i).Cost;
                    Structure(i).Cost=CostFunction(Offspring(2,:));
                end                
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            case 3 %Polynomial Mutation
                  Offspring=Structure(i).Decision;              
                  [proM,disM] = deal(1,12);
                  Site  = rand(1,nVar) < proM/nVar;
                  mu    = rand(1,nVar);
                  temp  = Site & mu<=0.5;
                  Offspring=Structure(i).Decision;
                  Offspring(temp) = Offspring(temp)+(VarMax(temp)-VarMin(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                                      (1-(Offspring(temp)-VarMin(temp))./(VarMax(temp)-VarMin(temp))).^(disM+1)).^(1/(disM+1))-1);
                  temp = Site & mu>0.5; 
                  Offspring(temp) = Offspring(temp)+(VarMax(temp)-VarMin(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                                      (1-(VarMax(temp)-Offspring(temp))./(VarMax(temp)-VarMin(temp))).^(disM+1)).^(1/(disM+1)));

                  Offspring = min(max(Offspring,VarMin),VarMax);  
                  Structure(i).Decision=Offspring;
                  Structure(i).PrevCost=Structure(i).Cost;
                  Structure(i).Cost=CostFunction(Offspring);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                   
            case 10 %Sin Cosin                
                a = 2;                                    
                costs=[Rep.Cost];
                cds=[Rep.CD];
                if all(cds==inf)
                    cds = ones(size(cds));
                else
                    cds(cds==inf) = 2*max(cds(cds~=inf));
                end
                s_costs=sum(costs,1);
                s_cds=sum(cds,1);
                P3=exp(-s_cds/sum(s_cds));
                P1=exp(-s_costs/sum(s_costs));
                P3=1-P3;
                P=Wcd*P3+Wsc*P1;
                P=P/sum(P);
                Selected_Sol=RouletteWheelSelection(P,1);

                mm=(Rep(Selected_Sol).Decision-Structure(i).Decision);

                r1=a-it*((a)/MaxIt); 

                r2=(2*pi)*randn();
                r3=2*randn;
                r4=randn();

                if r4<=0.5
                    newSol= Structure(i).Decision+r1*(sin(r2)*abs(r3*mm - ...
                        Structure(i).Decision));
                else                   
                    newSol= Structure(i).Decision+r1*(cos(r2)*(r3*mm - ...
                        Structure(i).Decision));
                end
             
                newSol = max(newSol, VarMin);
                newSol = min(newSol, VarMax);

                Structure(i).Decision=newSol;
                Structure(i).PrevCost=Structure(i).Cost;
                Structure(i).Cost=CostFunction(newSol);                       
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                
            case 2  % Extended line crossover  or BLX-0.25                
                costs=[Rep.Cost];
                cds=[Rep.CD];
                if all(cds==inf)
                    cds = ones(size(cds));
                else
                    cds(cds==inf) = 2*max(cds(cds~=inf));
                end
                s_costs=sum(costs,1);
                s_cds=sum(cds,1);
                P3=exp(-s_cds/sum(s_cds));
                P1=exp(-s_costs/sum(s_costs));
                P3=1-P3;
                P=(Wcd)*P3+(Wsc)*P1;
                P=P/sum(P);                
                Selected_Sol=RouletteWheelSelection(P,1);

                Parent1=Structure(i);
                Parent2=Rep(Selected_Sol(1));
                
                I=(Parent2.Decision-Parent1.Decision);
                r=unifrnd(-0.5,1.5, VarSize);

                Offspring=Parent1.Decision+r.*(Parent2.Decision-Parent1.Decision);
                         
                Offspring = min(max(Offspring,VarMin),VarMax);   
                
                Structure(i).Decision=Offspring;
                Structure(i).PrevCost=Structure(i).Cost;
                Structure(i).Cost=CostFunction(Offspring);               
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            case 20 % Modified GWO                
                costs=[Rep.Cost];
                cds=[Rep.CD];
                if all(cds==inf)
                    cds = ones(size(cds));
                else
                    cds(cds==inf) = 2*max(cds(cds~=inf));
                end
                s_costs=sum(costs,1);
                s_cds=sum(cds,1);
                P3=exp(-s_cds/sum(s_cds));
                P1=exp(-s_costs/sum(s_costs)); 
                P3=1-P3;
                P=Wcd*P3+Wsc*P1;
                P=P/sum(P);                
                Selected_Sol=RouletteWheelSelection(P,3);

                alpha=Rep(Selected_Sol(1));
                beta=Rep(Selected_Sol(1));
                delta=Rep(Selected_Sol(1));
                
                temp=zeros(1,nVar);
                a=2-it*((2)/MaxIt); 
                for kk=1:size(Structure(i).Decision,2)     

                    r1=unifrnd(-1,1); % 
                    r2=rand(); % 

                    A1=2*a*r1-a; % 
                    C1=2*r2; %
                    
                    D_alpha=abs(alpha.Decision(kk)-Structure(i).Decision(kk)); 
                    X1=(alpha.Decision(kk)-r1*D_alpha); 

                    r1=unifrnd(-1,1);
                    r2=rand();

                    A2=2*a*r1-a; 
                    C2=2*r2; % 

                    D_beta=abs(beta.Decision(kk)-Structure(i).Decision(kk)); 
                    X2=(beta.Decision(kk)-r1*D_beta);      

                    r1=unifrnd(-1,1);
                    r2=rand(); 

                    A3=2*a*r1-a; 
                    C3=2*r2; 

                    D_delta=abs(delta.Decision(kk)-Structure(i).Decision(kk)); 
                    X3=(delta.Decision(kk)-r1*D_delta);        
                    
                    temp(kk)=(X1+X2+X3)/3;             
                                        
                end
                Offspring=temp;
                Offspring = min(max(Offspring,VarMin),VarMax);   
                
                newC=CostFunction(Offspring);
                
                Structure(i).Decision=Offspring;
                Structure(i).PrevCost=Structure(i).Cost;
                Structure(i).Cost=CostFunction(Offspring);                     
                                             
                
        end  %END OF SWITCH
        
        Structure=DetermineDomination(Structure, i);
        Rep=MakeRepository(Structure, Rep);        
        Structure(i)=Update_Style(Structure(i), Selected_Style);
    end   
    
end  %End of function

