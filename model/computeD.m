function [Ds, concs] = computeD
close all;
clear

filename = 'gel_fitc_penetration_data.xlsx';
Data = xlsread(filename);
t = 60*45; % incubation time, in seconds 

Data = Data(1:201,:); % get only until the halfway point 


Ds = []; concs = {1,8};

Data_ctrl = Data(11:end,:);
[Dr,c] = fitData(Data_ctrl,3,t);
Ds = [Ds Dr];
concs{1,1} = c;

for i=4:5
    [Dr,c] = fitData(Data_ctrl, i, t);
    Ds = [Ds Dr];
    concs{1,i-2} = c;
end

[Dr,c] = fitData(Data_ctrl,6,t);
Ds = [Ds Dr];
concs{1,4} = c;

% This is for the ELAST curves: take away the initial part 
elast_index = 11; 
for i=7:10 
    Data_elast = Data(elast_index:end,:); 
    [Dr,c] = fitData(Data_elast, i, t);
    Ds = [Ds Dr];
    concs{1,i-2} = c; 
end


% Plotting 
norm_height = Data(:,2)/100; 

figure
hold on 
plot(norm_height(elast_index:end),fliplr(concs{1,1}*628.46/994),'b',norm_height,Data(:,3)/max(Data(:,3)),'b.')
plot(norm_height(elast_index:end),fliplr(concs{1,5}*1254/1272),'r',norm_height,Data(:,7)/max(Data(:,7)),'r.')
legend('Control (model)','Control (exp.)','ELAST (model)','ELAST (exp.)');
xlabel('Normalized depth')
ylabel('Normalized fluorescence')
title('Trial 1 (45 minutes)')

figure 
hold on
plot(norm_height(elast_index:end),fliplr(concs{1,2}*567.7/715.85),'b',norm_height,Data(:,4)/max(Data(:,4)),'b.')
plot(norm_height(elast_index:end),fliplr(concs{1,6}*1060.104/1087.785),'r',norm_height,Data(:,8)/max(Data(:,8)),'r.')
legend('Control (model)','Control (exp.)','ELAST (model)','ELAST (exp.)');
xlabel('Normalized depth')
ylabel('Normalized fluorescence')
title('Trial 2 (45 minutes)')

figure 
hold on
plot(norm_height(elast_index:end),fliplr(concs{1,3}*500.7/656.3),'b',norm_height,Data(:,5)/max(Data(:,5)),'b.')
plot(norm_height(elast_index:end),fliplr(concs{1,7}),'r',norm_height,Data(:,9)/max(Data(:,9)),'r.')
legend('Control (model)','Control (exp.)','ELAST (model)','ELAST (exp.)');
xlabel('Normalized depth')
ylabel('Normalized fluorescence')
title('Trial 3 (45 minutes)')

figure 
hold on
plot(norm_height(elast_index:end),fliplr(concs{1,4}*570/764),'b',norm_height,Data(:,6)/max(Data(:,6)),'b.')
plot(norm_height(elast_index:end),fliplr(concs{1,8}*1128.117/1142.569),'r',norm_height,Data(:,10)/max(Data(:,10)),'r.')
legend('Control (model)','Control (exp.)','ELAST (model)','ELAST (exp.)');
xlabel('Normalized depth')
ylabel('Normalized fluorescence')
title('Trial 4 (45 minutes)')

end


function [D1,c] = fitData(Data, column, real_t)
% Look here and change if want to fit all the data
index = length(Data(:,1))/1;
Depth_data = Data(:,1)*1000;
c_data_control_1 = Data(1:index,column)/max(Data(:,column)); % normalized
c_data_control_1 = fliplr(c_data_control_1');

% Input parameters of tissue
H = Depth_data(end); % in microns, denotes half the thickness of the tissue
D0 = 200; % initial guess, in microns^2/s

xdata = H - Depth_data; 
xdata = fliplr(xdata');
[D1, resnorm] = lsqcurvefit(@(D,xd) solvepde(D,xd,[H real_t]),D0,xdata,c_data_control_1);
D1
resnorm

% Solution profile (model) 
c = solvepde(D1, xdata, [H real_t]);
end