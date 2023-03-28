clear all
clc
close all

%% Constants and Set up
dim1 = 1e3;
dim2 = 1e3+1;

%parpool;

%% Launch and Arrival Porkchop
epoch = julian_date(2023, 0, 0, 0, 0, 0);

[VinfE, VinfM, launchT, transT] = CalcTransfer(epoch, 0, 1000, 0);

figure;
contourf(launchT, transT, VinfE.^2', 0:5:60)
caxis([0,60])
colormap('jet')
color = colorbar;
xlabel("Launch Wait (Days)")
ylabel("Transfer Time (Days)")
title("Earth Mars Transfer Departure C3")
ylabel(color, "C3 (km^2/s^2)")

figure;
contourf(launchT, transT, VinfM', 0:0.5:6)
caxis([0,6])
colormap('jet')
color = colorbar;
xlabel("Launch Wait (Days)")
ylabel("Transfer Time (Days)")
title("Earth Mars Transfer Arrival V_\infty")
ylabel(color, "V_\infty (km/s)")
