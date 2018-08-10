clear all
close all
load sim_40dots.mat

req_capacity=sum(debits);
BBU_capacity=20000;
BBU_cost=500000;
eq_used=zeros(nr_points,1);
pointcost=zeros(nr_points,1);

k=ceil(req_capacity/BBU_capacity);

nr_iterations=[1,2,5,10,20,50,100];
result=zeros(size(nr_iterations,2),1);

for it=1:size(nr_iterations,2)
    sim_costs=zeros(nr_iterations(it),1);

for j=1:nr_iterations(it)
    [idx,C,sumd,D] = kmeans(points,k);

    D=sqrt(D)/1000;

    total_cost=k*BBU_cost;

    for dots=1:nr_points
        [cost,eq_ref]=techtest(D(dots,idx(dots)),debits(dots));
        total_cost=total_cost+cost;
        eq_used(dots)=eq_ref;
        pointcost(dots)=cost;
    end
    sim_costs(j,1)=total_cost;
end
result(it)=min(sim_costs);
end

% figure
% hold on
% plot(points(:,1),points(:,2),'*');
% plot(C(:,1),C(:,2),'x'); 
