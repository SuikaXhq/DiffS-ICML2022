function subgroup = recover_full_index(subgroup, train_idx)

S = size(subgroup,2);
for s=1:S
    subgroup{s} = train_idx(subgroup{s});
end

end