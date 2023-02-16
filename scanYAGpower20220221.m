trial = ["V1", "V2", "V3", "V4", "V5"];
Ptables = {P100 P95 P90};
I500 = P100;
Itables = {I500 I1000 I1500 I2000};

%%
pedtime = 1:96;
abstime = 98 : size(P100,1);

%%

Ppdv = [];

for j = 1:size(Ptables,2)
    T = Ptables{j};
    for i = 1: size(trial,2)
        tri = trial(i);
        V = T.(tri);
        pdv = mean(V(pedtime));
        Ppdv(j, i) = pdv;
    end
end

Ipdv = [];

for j = 1:size(Itables,2)
    T = Itables{j};
    for i = 1: size(trial,2)
        tri = trial(i);
        V = T.(tri);
        pdv = mean(V(pedtime));
        Ipdv(j, i) = pdv;
    end
end







%%

figure(1)

for i = trial
    plot(P100.time, P100.(i))
    hold on;
end
title('YAG drive level : 100 %')
xlabel('time(s)')
ylabel('signal(V)')
legend(trial);
hold off;

figure(2)

for i = trial
    plot(P95.time, P95.(i))
    hold on;
end
title('YAG drive level : 95 %')
xlabel('time(s)')
ylabel('signal(V)')
legend(trial);
hold off;



figure(3)

for i = trial
    plot(P90.time, P90.(i))
    hold on;
end
title('YAG drive level : 90 %')
xlabel('time(s)')
ylabel('signal(V)')
legend(trial);
hold off;

%%

figure(4)
for i = 1 : size(Ptables,2)
    T = Ptables{i};
    for j = 1 : size(trial,2)
        tri = trial(j);
        V = T.(tri);
        pdv = Ppdv(i,j);
        plot(T.time, V/pdv);
        hold on;
    end
end
title('Normalized signals')
xlabel('time(s)')
ylabel('signal')
hold off;


%%

figure(5)

for i = trial
    T = I420;
    V = T.(i);
    plot(T.time, V)
    hold on;
end
title('I_abs = 420 µW')
xlabel('time(s)')
ylabel('signal(V)')
legend(trial);
hold off;

figure(6)

for i = trial
    T = I1000;
    V = T.(i);
    plot(T.time, V)
    hold on;
end
title('I_abs = 1000 µW')
xlabel('time(s)')
ylabel('signal(V)')
legend(trial);
hold off;

figure(7)

for i = trial
    T = I1500;
    V = T.(i);
    plot(T.time, V)
    hold on;
end
title('I_abs = 1500 µW')
xlabel('time(s)')
ylabel('signal(V)')
legend(trial);
hold off;

figure(8)

for i = trial
    T = I2000;
    V = T.(i);
    plot(T.time, V)
    hold on;
end
title('I_abs = 2000 µW')
xlabel('time(s)')
ylabel('signal(V)')
legend(trial);
hold off;

%%
figure(9)
for i = 1 : size(Itables,2)
    T = Itables{i};
    for j = 1 : size(trial,2)
        tri = trial(j);
        V = T.(tri);
        pdv = Ipdv(i,j);
        plot(T.time, V/pdv);
        hold on;
    end
end
title('Normalized signals')
xlabel('time(s)')
ylabel('signal')
hold off;