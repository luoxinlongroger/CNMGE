% 消耗时间对比图
NumOfEx = 1:1:34;
%large scale problem
CnmTime_CNMGE = [1.59E+01 1.60E+02 3.67E+02 6.83E+01 4.32E+00 5.78E+01 2.84E+01 1.08E+01 2.23E-01 1.34E+01 6.00E+00 4.96E+01 5.17E+01 2.03E+02 5.52E+00 1.87E+02 7.35E-01 2.15E+00 6.29E+01 5.89E+00 4.75E+01 1.44E+01 6.83E+00 3.98E+01 7.66E+01 3.92E+02 8.85E+01 3.58E+01 1.06E+01 2.04E+01 6.20E-01 7.20E+01 1.14E+02 9.11E+01];     
CnmTime_globalsearch = [18.1430 19.2532 3.7388 43.3965 3.4798 4.5220 4.3354 3.8359 3.2227 15.3921 5.1716 4.0257 3.8215 4.7188 3.8859 62.9534 10.3254 7.6031 4.5061 3.3021 3.3779 3.6557 3.3229 2.5044 5.4445 3.1651 22.1199 3.6483 3.5376 18.3366 11.6456 3.6256 2.7106 3.3033]; 
CnmTime_Couenne = [2.68E+04 1.77E+03 1.43E-01 2.42E+04 2.82E+01 2.55E+04 2.64E+04 4.14E-02 2.66E-03 3.26E-01 3.20E+00 2.13E-01 1.07E+02 9.05E-02 9.90E-02 2.62E+04 1.15E-01 9.68E-02 5.61E+02 2.58E+04 2.48E+00 4.10E+01 2.82E+04 3.49E-01 1.85E+01 1.23E+00 1.91E-01 2.67E+04 5.33E+01 7.34E-02 2.81E+04 1.60E+01 5.76E-02 1.12E+00];
CnmTime_CMA_ES = [2.98E+02 7.71E+01 1.51E+02 2.45E+02 1.84E+02 1.69E+02 1.84E+02 6.20E+01 5.91E+01 6.64E+02 1.61E+02 9.77E+01 1.77E+02 2.63E+02 5.75E+01 1.66E+02 2.83E+02 2.93E+02 2.66E+02 1.05E+02 6.64E+01 5.25E+01 9.50E+01 9.98E+01 2.44E+02 8.59E+01 1.22E+02 1.75E+02 2.52E+02 2.17E+02 1.65E+02 1.74E+02 4.47E+01 2.05E+02];
CnmTime_MCS = [2.01E+02 1.28E+04 1.10E+02 2.60E+02 1.25E+02 3.17E+02 5.45E+01 1.26E+04 1.26E+04 1.43E+04 1.44E+02 8.48E+02 2.18E+02 2.05E+02 1.27E+04 3.34E+02 1.22E+04 1.31E+04 1.59E+02 7.88E+01 6.79E+01 7.79E+01 9.22E+01 9.52E+01 7.05E+01 1.27E+04 1.41E+02 2.17E+02 7.52E+01 2.79E+02 1.25E+04 8.24E+01 1.01E+02 1.70E+02];
CnmTime_GLODS = [1.21E+02 1.14E+02 1.48E+02 1.35E+02 1.14E+02 1.27E+02 1.53E+02 1.13E+02 1.13E+02 1.30E+02 1.13E+02 1.19E+02 1.30E+02 1.22E+02 1.14E+02 8.91E+01 5.11E+01 5.16E+01 1.47E+02 1.30E+02 1.44E+02 1.53E+02 1.49E+02 6.67E+01 1.67E+02 1.36E+02 1.81E+02 0.01 7.09E+01 5.99E+01 0.01 1.45E+02 1.29E+02 1.30E+02];
CnmTime_VRBBO = [1.51E+01 4.81E+00 1.43E+01 1.36E+01 4.26E+00 3.82E+01 3.81E+00 3.59E+00 3.27E+00 1.71E+02 3.60E+00 3.98E+00 4.16E+00 1.95E+01 5.77E+00 3.93E+01 5.05E+00 5.08E+00 4.76E+00 3.36E+00 8.23E+00 3.69E+00 5.50E+00 4.05E+00 7.88E+00 3.38E+00 3.34E+00 7.51E+00 7.00E+00 1.18E+01 7.98E+00 3.43E+00 1.50E+01 3.98E+00];


figure(1)
plot(NumOfEx,CnmTime_CNMGE,'x-r');
%semilogy(NumOfEx,CnmTime_Ptc,'+-r');
hold on;
plot(NumOfEx,CnmTime_globalsearch,'d--b');
%semilogy(NumOfEx,CnmTime_Sqp,'d-b');
hold on;
plot(NumOfEx,CnmTime_Couenne,'^-k');
%semilogy(NumOfEx,CnmTime_BDF,'^-k');
hold on;
plot(NumOfEx,CnmTime_CMA_ES,'*-g');
% semilogy(NumOfEx,CnmTime_PM,'.-g');
hold on;
plot(NumOfEx,CnmTime_MCS,'.-c');
% semilogy(NumOfEx,CnmTime_PM,'.-c');
hold on

plot(NumOfEx,CnmTime_GLODS,'s-m');
% semilogy(NumOfEx,CnmTime_PM,'.-c');
hold on

plot(NumOfEx,CnmTime_VRBBO,'o--k');
% semilogy(NumOfEx,CnmTime_PM,'.-c');
hold off

%title('The error propagations of three algorithms','Color','k','FontSize',10);
% set(gca,'XTick',2:5:41);
% set(gca,'XTickLabel',{'2','3','5','6','7','8','9','11','12','13','15','16','17','18','19','20','21','22','23','24','25','27','28','30','32','33','34','35','36','37','38','39','40','41','42','44','45','46','47'});
set(gca,'YScale','log')     

xlabel('Problem','FontSize',10);
ylabel('Computational Time','FontSize',10);
legend('CNMGE','GlobalSearch','Couenne','CMA-ES','MCS','GLODS','VRBBO')

%small scale problem

NumOfEx = 35:1:68;
CnmTime_CNMGE = [1.98E+01 3.87E-01 1.40E+00 2.38E-01 1.97E+00 3.76E-02 4.42E-02 2.66E-02 3.58E-02 3.62E-02 5.67E-01 1.86E-01 1.00E-01 2.48E+00 3.93E-01 9.51E-02 9.84E-02 3.88E+00 3.82E-01 2.88E+00 8.16E-02 5.82E-02 1.83E-01 2.15E+00 3.14E+00 1.21E+00 6.77E-01 2.80E-01 5.12E-01 7.09E-02 2.13E+00 6.80E-02 5.32E-01 6.69E-02];      
CnmTime_globalsearch = [2.6565 2.8453 2.6163 3.7984 2.6378 5.2053 2.2066 2.4294 2.4533 2.1058 2.9943 4.1924 3.7953 4.1666 2.4326 2.3168 2.2858 3.1156 5.7896 2.7661 5.9145 2.3612 2.8148 2.9793 4.6011 3.3572 0.5485 3.5305 4.3224 2.7637 3.8752 2.4322 2.7524 2.7139 ]; 
CnmTime_Couenne = [1.14E-03 8.48E-03 2.80E+04 2.68E-03 2.55E-03 1.46E-03 8.57E-03 8.27E-01 1.05E-03 9.27E-01 7.47E-03 4.76E-03 2.39E-01 1.10E-03 2.43E-02 2.10E-01 9.38E-02 2.59E-02 2.10E-01 4.50E-01 7.88E-01 1.80E-02 7.61E-03 9.89E-01 1.49E-01 3.81E-01 2.31E-03 2.98E-02 3.62E-01 6.52E-03 7.46E-04 2.32E-03 6.78E-02 3.61E+00];
CnmTime_CMA_ES = [1.33 6.78E-01 5.70E-01 6.14E-01 3.71E-02 9.81E-01 6.35E-01 5.41E-01 5.49E-01 5.74E-01 5.37E-01 4.46 1.09 1.32 6.54E-01 5.50E-01 6.07E-01 7.34E-01 1.11 8.20E-01 8.23E-01 6.59E-01 6.26E-01 8.64E-01 9.21E-01 6.05E-01 5.67E-01 8.48E-01 5.97E-01 8.48E-01 5.69E-01 5.95E-01 4.10E-01 7.32E-01];
CnmTime_MCS = [6.38E-01 7.01E-01 3.26E-01 2.82E-01 6.13E-01 1.42E-01 4.91E-01 4.46E-01 2.48E-01 9.18E-02 2.30E-01 1.07 5.95E-01 6.62E-01 1.66E-01 1.15E-01 3.75E-01 5.65E-01 9.25E-01 3.93E-01 7.37E-01 2.57E-01 2.44E-01 7.22E-01 6.22E-01 5.86E-01 4.56E-01 4.70E-01 3.11E-01 1.59E-01 3.40E-01 4 2.64E-01 3.39E-01];
CnmTime_GLODS = [9.01E+00 5.80E+01 4.89E-01 5.82E+01 5.25E+01 5.48E+01 5.67E+01 5.36E+01 5.25E+01 5.23E+01 7.99E-01 2.24E+01 2.00E+01 5.37E+01 5.39E+01 5.37E+01 5.51E+01 5.32E+01 2.16E+01 3.23E+01 5.38E+01 5.39E+01 5.26E+01 5.33E+01 5.22E+01 1.92E+00 4.27E-01 5.26E+01 5.29E+01 3.28E+01 1.32E+00 5.43E+01 3.00E-01 1.86E+00];
CnmTime_VRBBO = [2.63E-01 3.74E-02 2.36E-02 4.92E-01 1.72E-01 2.91E-02 2.00E-02 1.79E-02 2.10E-02 3.77E-02 1.87E-02 3.11E-02 2.99E-02 2.03E-02 2.13E-02 2.26E-02 2.78E-02 1.74E-01 4.43E-02 3.76E-02 1.82E-01 1.84E-01 2.07E-02 2.44E-02 4.63E-02 2.17E-02 2.57E-02 2.35E-02 2.27E-02 3.01E-02 2.31E-02 2.21E-02 2.04E-02 2.69E-02];

figure(2)
plot(NumOfEx,CnmTime_CNMGE,'x-r');
%semilogy(NumOfEx,CnmTime_Ptc,'+-r');
hold on;
plot(NumOfEx,CnmTime_globalsearch,'d--b');
%semilogy(NumOfEx,CnmTime_Sqp,'d-b');
hold on;
plot(NumOfEx,CnmTime_Couenne,'^-k');
%semilogy(NumOfEx,CnmTime_BDF,'^-k');
hold on;
plot(NumOfEx,CnmTime_CMA_ES,'*-g');
% semilogy(NumOfEx,CnmTime_PM,'.-g');
hold on;
plot(NumOfEx,CnmTime_MCS,'.-c');
% semilogy(NumOfEx,CnmTime_PM,'.-c');
hold on

plot(NumOfEx,CnmTime_GLODS,'s-m');
% semilogy(NumOfEx,CnmTime_PM,'.-c');
hold on

plot(NumOfEx,CnmTime_VRBBO,'o--k');
% semilogy(NumOfEx,CnmTime_PM,'.-c');
hold off

%title('The error propagations of three algorithms','Color','k','FontSize',10);
% set(gca,'XTick',2:5:41);
% set(gca,'XTickLabel',{'2','3','5','6','7','8','9','11','12','13','15','16','17','18','19','20','21','22','23','24','25','27','28','30','32','33','34','35','36','37','38','39','40','41','42','44','45','46','47'});
set(gca,'YScale','log')     

xlabel('Problem','FontSize',10);
ylabel('Computational Time','FontSize',10);
legend('CNMGE','GlobalSearch','Couenne','CMA-ES','MCS','GLODS','VRBBO')


