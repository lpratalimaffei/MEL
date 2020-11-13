function [cell, body, new_species_number] = analyse_file(nominal_file, threshold)

if length(nominal_file.textdata) ==1
    A = strsplit(nominal_file.textdata{1,1});
else
    A = nominal_file.textdata;
end
% acquire species number
species_number = str2double(A{1,3});

% initialize some stuff
new_species_number = 0;
classifier = zeros(1,species_number);
counter = 0;

%initialize
body = [];

    for i = 1:species_number
        
        B = nominal_file.data(:,counter*3+2)>threshold;
        
        if sum(B) > 0
            
            body = [body, nominal_file.data(:,counter*3+1) ];
            body = [body, nominal_file.data(:,counter*3+2) ];
            body = [body, nominal_file.data(:,counter*3+3) ];
            
            new_species_number = new_species_number+1;
            classifier(1,i) = 1;
        else
            
        end
        
        counter = counter +1;
    end

%initialize
cell = {};    
cell{1,1} = A{1,1};
cell{1,2} = A{1,2};
cell{1,3} = num2str(new_species_number);

    for i = 1:species_number
       
        if classifier(1,i)==1
            cell{1,end+1}=A{1,3+i};
        end
    end
end