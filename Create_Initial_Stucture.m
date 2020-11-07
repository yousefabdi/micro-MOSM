
function Solutions=Create_Initial_Stucture()

    global ProblemSettings
    global SM_Settings
    
    CostFunction=ProblemSettings.CostFunction;
    VarSize=ProblemSettings.VarSize;
    nVar=ProblemSettings.nVar;
    VarMin=ProblemSettings.VarMin;
    VarMax=ProblemSettings.VarMax;
    nStyle=SM_Settings.nStyle;
    nPop=SM_Settings.nPop;
    
    candidate_solution.Decision=[];
    candidate_solution.Cost=[];
    candidate_solution.IsDominated=false;    
    candidate_solution.CD=0; 
    candidate_solution.Velocity=[];
      
    sol=repmat(candidate_solution,nPop,1);
    
    % Create Initial Candidate Solutions in sol Matrix
    for i=1:nPop

        sol(i).Decision(1,1:nVar)=unifrnd(VarMin,VarMax,VarSize);
                
        sol(i).Cost=CostFunction(sol(i).Decision);
        
        sol(i).Style=ones(1, nStyle);
        sol(i).PrevCost=sol(i).Cost;
        sol(i).CD=0;
                
    end
    Solutions = sol;
    
end % End of Function