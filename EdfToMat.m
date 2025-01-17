%EDF to mat conversion - https://www.mathworks.com/matlabcentral/fileexchange/38641-reading-and-saving-of-data-in-the-edf

[data, header]=ReadEDF('/MATLAB Drive/.edf');
data=cell2mat(data);

[data] = edfread('/MATLAB Drive/.edf');

save('.mat','header');
save('.mat','data');