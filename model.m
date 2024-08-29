clear;clc;close all;
H=10000/1000;beta=30*pi/180;Yc=H*tan(beta)/1000;
hold on
plot([Yc-3 Yc-1 Yc+1 Yc+3],[13.2 12.9 12.6 12.3],'ok', 'LineWidth', 10);
plot([Yc+3 Yc+1 Yc-1 Yc-3],[11.7 11.4 11.1 10.8],'ok', 'LineWidth', 10);
plot([Yc-12.5 Yc-10 Yc-7.5],[12.75 12.5 12.25 ],'ok', 'LineWidth', 10);
plot([Yc-7.5 Yc-10 Yc-12.5],[11.75 11.5 11.25],'ok', 'LineWidth', 10);
plot([Yc-12.5 Yc-12.5],[11.75 12.25],'ok', 'LineWidth', 10);
plot([Yc-10 Yc-5 Yc Yc+5],[12 12 12 12],'ok', 'LineWidth', 10);
xlim([Yc-20 Yc+20]);
ylim([10 15]);
xlabel('Range','FontSize',14);
ylabel('Azimuth','FontSize',14);
grid on
box off 

% 12.75 ,Yc-12.5
% 12.5,Yc-10
% 12.25,Yc-7.5
% 11.75,Yc-7.5
% 11.5,Yc-10
% 11.25,Yc-12.5
% 
% 13.2,Yc-3
% 12.9 ,Yc-1
% 12.6 ,Yc+1
% 12.3,Yc+3
% 11.7 ,Yc+3
% 11.4,Yc+1
% 11.1,Yc-1
% 10.8,Yc-3
% 
% 11.75 ,Yc-12.5
% 12.25,Yc-12.5
% 12,Yc-5
% 12 ,Yc
% 12 ,Yc-10
% 12 ,Yc+5

