function new_score = calculatePCs(data,coeff,center)
% Author: Jason Buenrostro, Stanford University
% calculates PC scores

% number of components
num_comp = size(coeff,2);

% center
if strcmp(center,'center')==1
    mean_data = nanmean(data);
    center_data = data - repmat(mean_data,[size(data,1) 1]);
else
    disp('Not centered!')
    center_data = data;
end

% calc score
clear new_score
for j = 1:num_comp
    % all data
    disp(j)
    for i = 1:size(data,2)
        new_score(i,j) = nansum(coeff(:,j).*center_data(:,i));
    end
end

end