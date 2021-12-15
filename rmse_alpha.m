function RMSE = rmse_alpha(v1, v2)

d = size(v1,2);
S = size(v1,1);
S_ = size(v2,1);
if S ~= S_
    RMSE = 0;
    return;
end

permlist = perms(1:S);
n = size(permlist,1);
RMSE = Inf;
for i=1:n
    RMSE_ = 1/sqrt(d)*norm(v1(permlist(i,:),:)-v2);
    if RMSE_ < RMSE
        RMSE = RMSE_;
    end
end

end