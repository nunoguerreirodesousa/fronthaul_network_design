clear all
l1=20;                                 %comprimento da area a considerar (km)
l2=20;                                 %largura da area a considerar (km)
density=1;                            %densidade de pontos por km^2

nr_points=ceil(density*(l1*l2));        %nr pontos a gerar tendo em conta densidade e area (arredonda para cima)

X=floor(rand(nr_points,1)*1000*l1);
Y=floor(rand(nr_points,1)*1000*l2);

possible_debits=[100,400,600,1000,1200,3000,6000,2000,];

points=[X,Y];

debits=zeros(nr_points,1);
for i=1:nr_points
    debits(i)=possible_debits(ceil(rand*size(possible_debits,2)));
end

x=[l1*500,l2*500];

figure
hold on
plot(points(:,1),points(:,2),'*');
plot(x(1,1),x(1,2),'x');

save sim_40dots.mat

%distance=dist2all(points,x);


function distance = dist2all(points,x)
% devolve array de distancia do ponto x a todos os pontos do array points
distance=zeros(size(points,1),1);
    
    for c= 1:size(points,1)
        distance(c,1)=sqrt(((points(c,1)-x(1,1))^2)+((points(c,2)-x(1,2))^2))/1000;
    end

end