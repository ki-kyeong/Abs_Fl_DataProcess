%% absorption curve

tarray = (1:2000);
figure(1)
plot(tarray, Data400.V1)
xlabel('time (µs)')
ylabel('signal (V)')
title('\Delta = 400 MHz')

%%

V200 = [Data200.V1 Data200.V2 Data200.V3 Data200.V4 Data200.V5];
V250 = [Data250.V1 Data250.V2 Data250.V3 Data250.V4 Data250.V5];
V300 = [Data300.V1 Data300.V2 Data300.V3 Data300.V4 Data300.V5];
V350 = [Data350.V1 Data350.V2 Data350.V3 Data350.V4 Data350.V5];
V400 = [Data400.V1 Data400.V2 Data400.V3 Data400.V4 Data400.V5];
V450 = [Data450.V1 Data450.V2 Data450.V3 Data450.V4 Data450.V5];
V500 = [Data500.V1 Data500.V2 Data500.V3 Data500.V4 Data500.V5];
V550 = [Data550.V1 Data550.V2 Data550.V3 Data550.V4 Data550.V5];
V600 = [Data600.V1 Data600.V2 Data600.V3 Data600.V4 Data600.V5];
V650 = [Data650.V1 Data650.V2 Data650.V3 Data650.V4 Data650.V5];
V700 = [Data700.V1 Data700.V2 Data700.V3 Data700.V4 Data700.V5];


%%
V200_1 = [];
V200_2 = [];
for i = 1:5
    V = V200(:,i);
    V = V(100:end);
    V200_1 = [V200_1 V];
    V200_2 = [V200_2 V/max(V)];
    S200_1 = sum(V200_1);
    S200_2 = sum(V200_2);
end

V250_1 = [];
V250_2 = [];
for i = 1:5
    V = V250(:,i);
    V = V(100:end);
    V250_1 = [V250_1 V];
    V250_2 = [V250_2 V/max(V)];
    S250_1 = sum(V250_1);
    S250_2 = sum(V250_2);
end

V300_1 = [];
V300_2 = [];
for i = 1:5
    V = V300(:,i);
    V = V(100:end);
    V300_1 = [V300_1 V];
    V300_2 = [V300_2 V/max(V)];
    S300_1 = sum(V300_1);
    S300_2 = sum(V300_2);
end

V350_1 = [];
V350_2 = [];
for i = 1:5
    V = V350(:,i);
    V = V(100:end);
    V350_1 = [V350_1 V];
    V350_2 = [V350_2 V/max(V)];
    S350_1 = sum(V350_1);
    S350_2 = sum(V350_2);
end

V400_1 = [];
V400_2 = [];
for i = 1:5
    V = V400(:,i);
    V = V(100:end);
    V400_1 = [V400_1 V];
    V400_2 = [V400_2 V/max(V)];
    S400_1 = sum(V400_1);
    S400_2 = sum(V400_2);
end

V450_1 = [];
V450_2 = [];
for i = 1:5
    V = V450(:,i);
    V = V(100:end);
    V450_1 = [V450_1 V];
    V450_2 = [V450_2 V/max(V)];
    S450_1 = sum(V450_1);
    S450_2 = sum(V450_2);
end

V500_1 = [];
V500_2 = [];
for i = 1:5
    V = V500(:,i);
    V = V(100:end);
    V500_1 = [V500_1 V];
    V500_2 = [V500_2 V/max(V)];
    S500_1 = sum(V500_1);
    S500_2 = sum(V500_2);
end

V550_1 = [];
V550_2 = [];
for i = 1:5
    V = V550(:,i);
    V = V(100:end);
    V550_1 = [V550_1 V];
    V550_2 = [V550_2 V/max(V)];
    S550_1 = sum(V550_1);
    S550_2 = sum(V550_2);
end

V600_1 = [];
V600_2 = [];
for i = 1:5
    V = V600(:,i);
    V = V(100:end);
    V600_1 = [V600_1 V];
    V600_2 = [V600_2 V/max(V)];
    S600_1 = sum(V600_1);
    S600_2 = sum(V600_2);
end

V650_1 = [];
V650_2 = [];
for i = 1:5
    V = V650(:,i);
    V = V(100:end);
    V650_1 = [V650_1 V];
    V650_2 = [V650_2 V/max(V)];
    S650_1 = sum(V650_1);
    S650_2 = sum(V650_2);
end

V700_1 = [];
V700_2 = [];
for i = 1:5
    V = V700(:,i);
    V = V(100:end);
    V700_1 = [V700_1 V];
    V700_2 = [V700_2 V/max(V)];
    S700_1 = sum(V700_1);
    S700_2 = sum(V700_2);
end

%%
freq = 200 : 50 : 700;

data1 = [mean(S200_1(2:end)) mean(S250_1) mean(S300_1) mean(S350_1) mean(S400_1) mean(S450_1) mean(S500_1) mean(S550_1) mean(S600_1) mean(S650_1) mean(S700_1)];
data2 = [mean(S200_2(2:end)) mean(S250_2) mean(S300_2) mean(S350_2) mean(S400_2) mean(S450_2) mean(S500_2) mean(S550_2) mean(S600_2) mean(S650_2) mean(S700_2)];

figure(2)
plot(freq, data1)
title('absorption w/o normalization')
xlabel('\Delta (MHz)')

figure(3)
plot(freq, data2)
title('absorption w/ normalization')
xlabel('\Delta (MHz)')

%%

figure(4)
h1 = waterfall(freq, tarray, [V200(:, 2) V250(:, 2) V300(:, 2) V350(:, 2) ...
    V400(:, 2) V450(:, 2) V500(:, 2) V550(:, 2) V650(:, 2) V700(:, 2) ] );

hold on;
