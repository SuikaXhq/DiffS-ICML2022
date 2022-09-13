S = size(subgroup,2);
M = size(division_index,1);
color = zeros(1,M);
for s=1:S
    for m=subgroup{s}
        color(m) = s;
    end
end
color(48)=[];
file = fopen('real_color_kmns.txt','a');
fprintf(file, sprintf('# k=%d\ncolor_exp = ', K));
for m=1:M-1
    fprintf(file, sprintf('%d,', color(m)));
end
fprintf(file, sprintf('%d\n', color(m)));