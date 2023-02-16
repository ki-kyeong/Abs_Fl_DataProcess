path = "2022-12-27_Vac.csv";
Data = readtable(path);
figure(1);
x=cell2mat(Data{:,1});
y = Data{:,2};

figure(1)
plot( y)
