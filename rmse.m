function RMSE = rmse(v1, v2, test_idx)
if nargin < 3
    test_idx = 1:size(v1,1);
end

d = size(v1,2);
RMSE = 1/sqrt(d)*norm(v1(test_idx,:)-v2(test_idx,:));

end