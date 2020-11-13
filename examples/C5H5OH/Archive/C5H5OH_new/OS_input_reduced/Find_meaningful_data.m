%% FIND MEANINGFUL SPECIES TO LOOK AT %%
clear; clc;
% Import files
List_of_Exp = importdata('Path_to_Exp_Datasets.txt');
List_of_Inp = importdata('Path_to_OS_inputs.txt');
% Define here the quantities of interest
T  = [600 800 1000 1200 1400 1600 1800 2000 2200 2400];
P  = [0.01 0.1 1 10 100];
SP = {'P1', 'P2_L', 'P3', 'P4', 'R', 'W3', 'W4_L'};


% initialize
cell = {};

for i = 1:size(List_of_Exp,1)
    
    A = strsplit(List_of_Exp{i,1},'/');
    
    single_exp = importdata(List_of_Exp{i,1});
    
    % collect time
    Time_temp  = [single_exp.data(1,1), single_exp.data(end,1), (single_exp.data(end,1)-single_exp.data(1,1))/2];
    
    % choose species
    if strcmp(A{1,1},SP{1,1})
            col = 1;
    elseif strcmp(A{1,1},SP{1,2})
            col = 2;
    elseif strcmp(A{1,1},SP{1,3})
            col = 3;
    elseif strcmp(A{1,1},SP{1,4})   
            col = 4;
    elseif strcmp(A{1,1},SP{1,5})   
            col = 5;
    elseif strcmp(A{1,1},SP{1,6})   
            col = 6;
    elseif strcmp(A{1,1},SP{1,7})   
            col = 7;
    end
    
    % choose temperature
    if strcmp(A{1,3},[num2str(T(1)) 'K'])
            row = 1;
    elseif strcmp(A{1,3},[num2str(T(2)) 'K'])
            row = 2;
    elseif strcmp(A{1,3},[num2str(T(3)) 'K'])
            row = 3;
    elseif strcmp(A{1,3},[num2str(T(4)) 'K'])
            row = 4;
    elseif strcmp(A{1,3},[num2str(T(5)) 'K'])
            row = 5;
    elseif strcmp(A{1,3},[num2str(T(6)) 'K'])
            row = 6;
    elseif strcmp(A{1,3},[num2str(T(7)) 'K'])
            row = 7;
    elseif strcmp(A{1,3},[num2str(T(8)) 'K'])
            row = 8;
    elseif strcmp(A{1,3},[num2str(T(9)) 'K'])
            row = 9;
    elseif strcmp(A{1,3},[num2str(T(10)) 'K'])
            row = 10;
    end
    
    % choose pressure
    if strcmp(A{1,2},[num2str(P(1),'%10.2f') 'atm'])
            line = 1;
    elseif strcmp(A{1,2},[num2str(P(2),'%10.1f') 'atm'])
            line = 2;
    elseif strcmp(A{1,2},[num2str(P(3),'%10.1f') 'atm'])
            line = 3;
    elseif strcmp(A{1,2},[num2str(P(4),'%10.1f') 'atm'])
            line = 4;
    elseif strcmp(A{1,2},[num2str(P(5),'%10.1f') 'atm'])
            line = 5;
    end
    
    cell{row,col}(line,:) = Time_temp;
end

% % ciao
% scatter(P, cell{1,1}(:,3), 'dk','filled')
% % ciaone
% CHECK = errorbar(1000 ./ Cond_P(:,14),Cond_P(:,17),Cond_P(:,17)*0.5,Cond_P(:,17)*0.5, color_p{p},'LineStyle','none');
% set(get(get(CHECK(),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

%%
color_marker = {'dk', 'dr', 'db', 'dc', 'dm', 'dg', 'dy'};

for i = 1:size(cell,1)
    
    figure()
    title([num2str(T(i)) ' K'])
    
    set(gca, 'Xscale', 'log');
    set(gca, 'Yscale', 'log');
    
    xlabel('Pressure [atm]') 
    ylabel('Average reactivity time [s]') 
    set(gca,'FontSize',20)
    
    hold on;
    
    for j=1:size(cell,2)
        
        scatter(P, cell{i,j}(:,3), color_marker{j},'filled'); set(gca, 'Yscale', 'log');
        CHECK = errorbar(P, cell{i,j}(:,3),cell{i,j}(:,1),cell{i,j}(:,2), 'k','LineStyle','none');
        set(get(get(CHECK(),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    end
    
    lgd = legend(SP , 'FontSize', 12, 'Location' , 'Best');
    
end