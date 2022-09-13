function [subgroup_est, theta] = estimate_groups(subgroup_init, alpha, theta)

subgroup_est = subgroup_init;
M = size(theta, 1);
S = size(alpha, 1);
isinitialized = zeros(1,M);
for s=1:S
    for member=subgroup_init{s}
        isinitialized(member) = 1;
    end
end

for i=1:M
    if ~isinitialized(i)
        diff = alpha - theta(i,:);
        min_dist = Inf;
        for s=1:S
            if norm(diff(s,:)) < min_dist
                min_dist = norm(diff(s,:));
                subgroup_idx_est = s;
            end
        end
        subgroup_est{subgroup_idx_est}(end+1) = i;
    end
end

for s=1:S
    for member=subgroup_est{s}
        theta(member, :) = alpha(s,:);
    end
end

end