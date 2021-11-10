ncolor = 20;
cocV1 = cbrewer('seq', 'YlGnBu', ncolor);

a = linspace(1,500,1000);

figure
for i=1:ncolor
    v = ones(1000,1)*i;
    plot(a,v, 'color', cocV1(i,:))
    hold all
end

%%

ncolor = 5;

colPl = cbrewer('seq', 'RdPu',30);
cocNaiveV1 = colPl(13:17,:);
cocEarlyV1 = colPl(19:24,:);
cocAdultV1 = colPl(26:30,:);

cocV1 = [cocNaiveV1; cocEarlyV1; cocAdultV1];
a = linspace(1,500,1000);

figure
for i=1:ncolor*3
    v = ones(1000,1)*i;
    plot(a,v, 'color', cocV1(i,:))
    hold all
end