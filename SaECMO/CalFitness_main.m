function Fitness = CalFitness_main(Population,VAR)
%--------------------------------------------------------------------------
PopObj = Population.objs;
PopCon = Population.cons;
N = length(Population);
CV = sum(max(0,PopCon),2);
CV_sorted = sort(CV);
t = N/4;%threhold
M_VAR = max(N-floor(VAR * N), N/4);

CV_m = CV_sorted(M_VAR);
%% Detect the dominance relation between each two solutions
Dominate = false(N);
for i = 1 : N-1
    k_cv_i = logical(CV(i)<=CV_m);
    for j = i+1 : N
        k_obj = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
        k_cv_j = logical(CV(j)<=CV_m);             
        if k_cv_i - k_cv_j == 0        
            if k_cv_i == 1
                if k_obj == 1
                    Dominate(i,j) = true;
                elseif k_obj == -1
                    Dominate(j,i) = true;
                end
            else
                if CV(i)<CV(j)
                    Dominate(i,j) = true;
                end
            end            
        elseif k_cv_i ==1
            Dominate(i,j) = true;
        else
            Dominate(j,i) = true;          
        end       
    end
end

%% Calculate S(i)
S = sum(Dominate,2);

%% Calculate R(i)
R = zeros(1,N);
for i = 1 : N
    R(i) = sum(S(Dominate(:,i)));
end

%% Calculate D(i)
Distance = pdist2(PopObj,PopObj);
Distance(logical(eye(length(Distance)))) = inf;
Distance = sort(Distance,2); %每一行从小到大排列
D = 1./(Distance(:,floor(sqrt(N)))+2);


%% Calculate the fitnesses
Fitness = R + D';
end