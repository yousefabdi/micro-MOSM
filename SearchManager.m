%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Matlab Source Code for 
%%%%%%%%%%%  Configurations of Micro Multi-objective Search Manager
%%%%%%%%%%%  Programmed By: Yousef Abdi
%%%%%%%%%%%  E-mail: yousef.abdi@gmail.com, y.abdi@tabrizu.ac.ir
%%%%%%%%%%%  November 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
format longE;

for Problem=[1 2 3 4 6]
    %%  Globalization of Parameters and Settings
    prob=num2str(Problem);
    CostFunction=str2func(['@(x) ZDT' prob '(x)']);  %Cost Function
    load(['Problem_PFs\ZDT' prob]);
    problemName=['ZDT' prob ];
    
    VarMin=[]; VarMax=[];
    
    switch Problem
        case num2cell([1 2 3])
            nVar=30;            
            VarMin(1:nVar)=0;
            VarMax(1:nVar)=1;
        case 4
            nVar=10;
            VarMin(1:nVar)=0;
            VarMax(1:nVar)=1;
           VarMin(2:nVar)=-5;
           VarMax(2:nVar)=5;
        case 6
            nVar=10;
            VarMin(1:nVar)=0;
            VarMax(1:nVar)=1;          
    end
    
    global ProblemSettings;
    global SM_Settings;

    VarSize=[1 nVar];   % Decision Variables Matrix Size

    %% SEARCH MANAGER Parameters

    MaxIt=2000;         % Maximum Number of Iterations
    TotalRun=1;
    nPop=5;            % Population Size
    
    % Number of Optimization Styles that want to be hybridized
    % Note that the nStyle variable should be set to 3 for micro-MOSM-C2
    nStyle=2;           
    
    %% Initialization
    ProblemSettings.CostFunction=CostFunction;
    ProblemSettings.nVar=nVar;
    ProblemSettings.VarSize=VarSize;
    ProblemSettings.VarMin=VarMin;
    ProblemSettings.VarMax=VarMax;
    ProblemSettings.MaxIt=MaxIt;

    SM_Settings.MaxIt=MaxIt;
    SM_Settings.nPop=nPop;
    SM_Settings.nRep=100;
    SM_Settings.nStyle=nStyle;

    IGDs=zeros(MaxIt, TotalRun);
    HVs=zeros(MaxIt, TotalRun);

    for run=1:TotalRun

        it=0;
        Structure=Create_Initial_Stucture(); 

        %% Determine Domination & Create Initial Repository
        Rep = [];
        Structure=DetermineDomination(Structure, 0);
        Rep = MakeRepository(Structure, Rep);

        %% Search Manager Main Loop
        for it=1:MaxIt    

            Prev_Rep=[Rep.Cost];

            % Movement Operation
            Structure=Movement(Structure, Rep, it);        

            % Dominations ...
            Structure=DetermineDomination(Structure, 0);

            Rep=MakeRepository(Structure, Rep);        
     
            rep_cos=[Rep.Cost];
            
            IGDs(it,run)=IGD([Rep.Cost]',PF);
            HVs(it,run)=HV([Rep.Cost]',PF);
            disp(['Problem= ' num2str(Problem) ' : Run= ' num2str(run) ' : Iteraion ' num2str(it) ': Number of Repository Members = ' num2str(numel(Rep)) ': HV = ' num2str(HVs(it,run))]);

            %Plot Costs
            figure(1);
            PlotCosts2([Rep.Cost]);
            pause(0.0001);          
       end
     
    end

end
    

