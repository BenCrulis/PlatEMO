function [Population,FrontNo,Fitness,D] = EnvironmentalSelection(Population,N,IdealPoint)
% The environmental selection of C-CLOP


% This function is written by Ben Crulis

    %% let's round the objective values
    %objs = round(Population.objs, 4);
    objs = Population.objs;
    
    %% Non-dominated sorting
    [FrontNo, Last] = NDSort(objs,Population.cons,N);
    %[FrontNo, Last] = NDSort(objs, N);
    Front = FrontNo == 1;

    [nInd, ~]    = size(objs);
    
    Infeas = any(Population.cons > 0, 2);
    ConsF = (Population.cons > 0) .* Population.cons;
    Cons = sum(ConsF, 2);
        
    %% Compute the normalized scores
    F1 = objs(Front,:);
    
    Nadir = max(F1, [], 1);
    %Nadir = mean(F1, 1); % using max is unstable for some problems such as IDTLZ1
    
    Normalization = Nadir - IdealPoint;
    Normalization(Normalization == 0) = 1;

    F = (objs - repmat(IdealPoint, nInd,1))./repmat(Normalization, nInd,1);
        
    %% compute the fitness as the distance to the ideal point

    Fitness = squeeze(sum(F.^2,2));
    %Fitness = sum(F, 2);
    
    %% compute the distance between solutions
    
    %F = WS(F + 1e-7); % transform the vectors using MOEA/D-AWA WS transformation
    
    D = pdist2(F,F,'cosine');
    
    D(isnan(D)) = 0;
    
    D(logical(eye(nInd))) = inf;
    
    %% Truncate solutions
    Keep = false(1,nInd);
    Keep(FrontNo <= Last) = true;
    
    %Keep(Infeas) = false;
    
    %if sum(Keep) < N
        %Keep(~Infeas) = true;
    %end
    
    Corner = FindCornerSolutions(F);
    %Corner = [];
    
    while sum(Keep) > N
        Remain   = find(Keep);
        [B, I]   = min(D(Remain,Remain), [],2);
        [~,Second] = min(B(:,1));
        First = Remain(I(Second,1));
        Second = Remain(Second);

        ToDel = First;

        FirstIsCorner = ismember(First, Corner);
        SecondIsCorner = ismember(Second, Corner);
        
        consFirst = Cons(First);
        consSecond = Cons(Second);
        
        if ((consFirst > 0) || (consSecond > 0)) && (consFirst ~= consSecond)
            if consFirst < consSecond
                ToDel = Second;
            else
                ToDel = First;
            end    
        elseif FirstIsCorner && ~SecondIsCorner
            ToDel = Second;
        elseif ~FirstIsCorner && SecondIsCorner
            ToDel = First;
        elseif Fitness(First) < Fitness(Second)
            ToDel = Second;
        end
        Keep(ToDel) = false;
    end
        
    %% Population for next generation
    Population = Population(Keep);
    FrontNo    = FrontNo(Keep);
    %Fitness = sum(mink(D, size(objs,2), 2), 2); % use angle fitness
    Fitness    = Fitness(Keep);
    D = D(Keep,Keep);
end