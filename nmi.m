function [NMI, is_perfect] = nmi(subgroup, subgroup_est)
S = size(subgroup,2);
S_est = size(subgroup_est,2);
for s=1:S
    if size(subgroup{s},1)==0
        subgroup(s) = [];
    end
end
for s=1:S_est
    if size(subgroup_est{s},1)==0
        subgroup_est(s) = [];
    end
end

S = size(subgroup,2);
S_est = size(subgroup_est,2);
M_s = zeros(1,S);
M_s_est = zeros(1,S_est);
for s=1:S
    M_s(s) = size(subgroup{s}, 2);
    subgroup{s} = sort(subgroup{s});
end
for s=1:S_est
    M_s_est(s) = size(subgroup_est{s}, 2);
    subgroup_est{s} = sort(subgroup_est{s});
end
M = sum(M_s);

I = 0;
for s=1:S
    for s_=1:S_est
        K = size(intersect(subgroup{s},subgroup_est{s_}),2);
        if K~=0
            I = I + K/M*log(M*K/M_s(s)/M_s_est(s_));
        end
    end
end

H = -sum(M_s/M.*log(M_s/M));
H_est = -sum(M_s_est/M.*log(M_s_est/M));

NMI = 2*I/(H+H_est);

%% Is perfect?
is_perfect = true;
if S == S_est
    for s=1:S
        flag = false;
        for s_=1:S
            if M_s(s) == M_s_est(s_)
                if sum(subgroup{s}==subgroup_est{s_})==M_s(s)
                    flag = true;
                    break;
                end
            end
        end
        if ~flag
            is_perfect = false;
            break;
        end
    end
else
    is_perfect = false;
end

end