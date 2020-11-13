%% ciao %%
clear; clc;

List_of_Exp = importdata('Path_to_Exp_Datasets.txt');
List_of_Inp = importdata('Path_to_OS_inputs.txt');

tot_profiles = 0;

for i = 1:size(List_of_Exp,1)
    
    if i ~= 1       
       clear single_exp;
    end
    
    
    % create directories
    create_the_directories(List_of_Exp{i,1});
    % copy OS++ input file
    unix(['cp ' List_of_Inp{i,1} ' OUTPUT/' List_of_Inp{i,1}  ]);
    
    single_exp = importdata(List_of_Exp{i,1});
    [cell_header, body, new_species_number] = analyse_file(single_exp, 1.0e-9);
    
    tot_profiles = tot_profiles+new_species_number;
    %create file name
    filename = fullfile(['OUTPUT/' List_of_Exp{i,1}]);
    %write header
    fid = fopen(filename, 'wt+');
        fprintf(fid, '%s\t%s\t%s\n', strjoin(cell_header,'\t')); % header
        fprintf(fid, '\n');
    fclose(fid);
    %write body
    dlmwrite(filename,body,'delimiter','\t','precision',['%10.',num2str(12),'e'],'-append');
end

fprintf(['The total number of retained profile is ' num2str(tot_profiles) '\n']);
fprintf(['The removed ones are ' num2str(100-tot_profiles*100/(size(List_of_Exp,1)*4)) ' percent \n' ]);