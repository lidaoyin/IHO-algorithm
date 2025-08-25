function[Best_score,Best_pos,IHO_curve]=IHO(lowerbound,upperbound,fitness,SearchAgents,dimension,Max_iterations)



lowerbound=ones(1,dimension).*(lowerbound);                              % Lower limit for variables
upperbound=ones(1,dimension).*(upperbound);                              % Upper limit for variables

%% Initialization
for i=1:dimension
    X(:,i) = lowerbound(i)+rand(SearchAgents,1).*(upperbound(i) - lowerbound(i));                          % Initial population
end

for i =1:SearchAgents
    L=X(i,:);
    fit(i)=fitness(L);
end


%% Main Loop
for t=1:Max_iterations
    %% Update the Best Condidate Solution
    [best , location]=min(fit);
    if t==1
        Xbest=X(location,:);                                  % Optimal location
        fbest=best;                                           % The optimization objective function
    elseif best<fbest
        fbest=best;
        Xbest=X(location,:);
    end
    
    Dominant_hippopotamus=Xbest;
    
    %population initialization
    [fit_sort,fit_sort_index]=sort(fit);
    X_sort=X(fit_sort_index,:);
    
    
    for i=1:floor(0.3*SearchAgents)
        
        fit_esc(i)=fit_sort(i);
        X_esc(i,:)=X_sort(i,:);

        f_esc_before(t,i)=fit_esc(i);
        
        LO_LOCAL=(lowerbound./t);
        HI_LOCAL=(upperbound./t);
        Alfa{1,:}= 2*rand(1,dimension)-1;
        Alfa{2,:}= rand(1,1);
        Alfa{3,:}=randn;
        D=Alfa{randi([1,3],1,1),:};
        X_P4(i,:)=X_esc(i,:)+(rand(1,1)).*(LO_LOCAL+D.* (HI_LOCAL-LO_LOCAL));
        X_P4(i,:) = min(max(X_P4(i,:),lowerbound),upperbound);

        L=X_P4(i,:);
        F_P4(i)=fitness(L);
        if(F_P4(i)<fit_esc(i))
            X_esc(i,:) = X_P4(i,:);
            fit_esc(i) = F_P4(i);
        end

        fit_sort(i)=fit_esc(i);
        X_sort(i,:)=X_esc(i,:);
        
        f_esc_after(t,i)=fit_esc(i);
        f_esc_rec(t,i)=f_esc_before(t,i)-f_esc_after(t,i);        
        
    end
    
    
    [fbest,location]=min(fit_sort);
    Xbest=X_sort(location,:);
    Dominant_hippopotamus=Xbest;
    
    %Defense against predators
    for i=floor(0.3*SearchAgents)+1:floor(0.6*SearchAgents)
        
        R1=0.15*(rand+1);
        R2=0.4*randi([1,2],1,1);
        Alfa{1,:}=R1*(2*rand(1,dimension)-1);%near position
        Alfa{2,:}=R2*(2*rand(1,dimension)-1);%remote location
        D1=randi([1,2],1,1);%Randomly selecting a close location or a distant location
        D_P=Alfa{D1,:};
        
        predator=X_sort(i,:)+D_P.*(upperbound-lowerbound); 
        predator=min(max(predator,lowerbound),upperbound);
        F_predator=fitness(predator);
        
        f_p1_before(t,i)=fit_sort(i);
        
        if fit_sort(i)> F_predator
            fit_sort(i)=F_predator;
            X_sort(i,:)=predator;
        else
            if D1==1
                x_predator1_new=X_sort(i,:)+rand*(Dominant_hippopotamus-abs(X_sort(i,:))) - rand*(predator-abs(X_sort(i,:)));
            else
                x_predator1_new=X_sort(i,:)-rand(1,dimension).*(predator-X_sort(i,:));
            end
            x_predator1_new= min(max(x_predator1_new,lowerbound),upperbound);
            F_x_new=fitness(x_predator1_new);
            if(F_x_new<fit_sort(i))
                fit_sort(i)=F_x_new;
                X_sort(i,:)=x_predator1_new;
            end
        end
        f_p1_after(t,i)=fit_sort(i);
        f_p1_rec(t,i)=f_p1_before(t,i)-f_p1_after(t,i);
    end
    
    
    [fbest,location]=min(fit_sort);
    Xbest=X_sort(location,:);
    Dominant_hippopotamus=Xbest;
    
    for i=floor(0.6*SearchAgents)+1:SearchAgents
    
        f_1_before(t,i)=fit_sort(i);
        
        I1=randi([1,2],1,1);
        I2=randi([1,2],1,1);
        Ip1=randi([0,1],1,2);
        RandGroupNumber=randperm(SearchAgents,1);
        RandGroup=randperm(SearchAgents,RandGroupNumber);

        % Mean of Random Group
        MeanGroup=mean(X_sort(RandGroup,:)).*(length(RandGroup)~=1)+X_sort(RandGroup(1,1),:)*(length(RandGroup)==1);
        Alfa{1,:}=(I2*rand(1,dimension)+(~Ip1(1)));
        Alfa{2,:}= 2*rand(1,dimension)-1;
        Alfa{3,:}= rand(1,dimension);
        Alfa{4,:}= (I1*rand(1,dimension)+(~Ip1(2)));
        Alfa{5,:}=rand;
        A=Alfa{randi([1,5],1,1),:};
        B=Alfa{randi([1,5],1,1),:};
        X_P1(i,:)=X_sort(i,:)+rand(1,1).*(Dominant_hippopotamus-I1.*X_sort(i,:));
        T=exp(-t/Max_iterations);
        if T>0.6
            X_P2(i,:)=X_sort(i,:)+A.*(Dominant_hippopotamus-I2.*MeanGroup);
        else
            if rand()>0.5
                X_P2(i,:)=X_sort(i,:)+B.*(MeanGroup-Dominant_hippopotamus);
            else
                X_P2(i,:)=((upperbound-lowerbound)*rand+lowerbound);
            end

        end
        X_P2(i,:) = min(max(X_P2(i,:),lowerbound),upperbound);

        L=X_P1(i,:);
        F_P1(i)=fitness(L);
        if(F_P1(i)<fit_sort(i))
            X_sort(i,:) = X_P1(i,:);
            fit_sort(i) = F_P1(i);
        end

        L2=X_P2(i,:);
        F_P2(i)=fitness(L2);
        if(F_P2(i)<fit_sort(i))
            X_sort(i,:) = X_P2(i,:);
            fit_sort(i) = F_P2(i);
        end


        f_1_after(t,i)=fit_sort(i);
        f_1_rec(t,i)=f_1_before(t,i)-f_1_after(t,i);
        
    end
    
    
    X=X_sort;
    fit=fit_sort;
    
    fbest=min(fit);
    best_so_far(t)=fbest;
%     disp(['Iteration ' num2str(t) ': Best Cost = ' num2str(best_so_far(t))]);




end
Best_score=fbest;
Best_pos=Xbest;
IHO_curve=best_so_far;
    

end