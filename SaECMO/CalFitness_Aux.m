function Fitness = CalFitness_Aux(Population,VAR)
% Calculate the fitness of each solution

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
PopObj = Population.objs;
PopCon = Population.cons;
N = length(Population);
CV = sum(max(0,PopCon),2);
CV_sorted = sort(CV);
ifindex = find(CV_sorted>0);
if_N = length(ifindex);
c = round((1-VAR) * if_N);
a = max(c, 1);
if if_N < a
    CV_m = 0;
else
    M_VAR = ifindex(a);
    M_VAR = max(M_VAR, 1);
    CV_m = CV_sorted(M_VAR);   
end

z = [PopObj,CV];
%% Detect the dominance relation between each two solutions
Dominate = false(N);
for i = 1 : N-1
    k_cv_i = logical(CV(i)<=CV_m);
    for j = i+1 : N
        k_obj = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
        k_cv_j = logical(CV(j)<=CV_m);        
        k = any(z(i,:)<z(j,:)) - any(z(i,:)>z(j,:));
        if k_cv_i - k_cv_j == 0
                if k_obj == 1
                    Dominate(i,j) = true;
                elseif k_obj == -1
                    Dominate(j,i) = true;
                end      
    elseif k_cv_i ==1
        if k_obj == 1         
            Dominate(i,j) = true;
        end
        elseif k_obj == -1
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