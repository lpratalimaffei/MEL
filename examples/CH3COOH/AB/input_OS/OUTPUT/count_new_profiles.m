%% ciao %%

addpath('/Users/andrea/Dropbox/LPM/CH3COOH_new/input_OS/')
clear; clc;

List_of_Exp = importdata('Path_to_Exp_Datasets_new.txt');
List_of_Inp = importdata('Path_to_OS_inputs_new.txt');

tot_profiles = 0;

for i = 1:size(List_of_Exp,1)
    
    if i ~= 1       
       clear single_exp;
    end
    
    
    single_exp = importdata(List_of_Exp{i,1});
    
    if length(single_exp.textdata) ==1
        A = strsplit(single_exp.textdata{1,1});
    else
        A = single_exp.textdata;
    end
    
    new_species_number = str2double(A{1,3});
    
    tot_profiles = tot_profiles+new_species_number;
    
end

fprintf(['The total number of retained profile is ' num2str(tot_profiles) '\n']);
fprintf(['The removed ones are ' num2str(100-tot_profiles*100/(208*4)) ' percent of the total initial \n' ]);