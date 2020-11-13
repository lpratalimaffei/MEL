function create_the_directories(string)

A = strsplit(string,'/');

dir1 = A{1,1};
dir2 = [A{1,1} '/' A{1,2}];
dir3 = [A{1,1} '/' A{1,2} '/' A{1,3}];

if ~exist(['OUTPUT/' dir1], 'dir')
       mkdir(['OUTPUT/' dir1])
end

if ~exist(['OUTPUT/' dir2], 'dir')
       mkdir(['OUTPUT/' dir2])
end

if ~exist(['OUTPUT/' dir3], 'dir')
       mkdir(['OUTPUT/' dir3])
end


end