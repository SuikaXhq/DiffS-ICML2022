temp = textread('real_data/climdiv-tmpcdv-v1.0.0-20200804');
phdi = textread('real_data/climdiv-phdidv-v1.0.0-20200804');
pcpn = textread('real_data/climdiv-pcpndv-v1.0.0-20200804');
pdsi = textread('real_data/climdiv-pdsidv-v1.0.0-20200804');
zndx = textread('real_data/climdiv-zndxdv-v1.0.0-20200804');

division_index_full = floor(temp(:,1)/1000000);

temp = temp(1:43344, 2:end);
phdi = phdi(1:43344, 2:end);
pcpn = pcpn(1:43344, 2:end);
pdsi = pdsi(1:43344, 2:end);
zndx = zndx(1:43344, 2:end);

% count total number of divisions
division_index = unique(division_index_full);
division_index = division_index(1:344);
M = size(division_index,1);
n = sum(division_index_full'==division_index, 2);
N = sum(n);

% beta = (summer, fall, winter, PDSI, PHDI)
% theta = (Intercept(constant 1), PCPN, ZNDX)
Y = cell(1,M);
X = cell(1,M);
Z = cell(1,M);
for i=1:M
    Y{i} = zeros(n(i)*12, 1);
    X{i} = zeros(n(i)*12, 5);
    Z{i} = zeros(n(i)*12, 3);
    for k=1:n(i)
        for m=1:12
            offset =12*(k-1) + m;
            index = sum(n(1:i-1)) + k;
            Y{i}(offset) = temp(index, m);
            if m==12 || m==1 || m==2
                X{i}(offset, 3) = 1;
            elseif m>=6 && m<=8
                X{i}(offset, 1) = 1;
            elseif m>=9 && m<=11
                X{i}(offset, 2) = 1;
            end
            X{i}(offset, 4) = pdsi(index, m);
            X{i}(offset, 5) = phdi(index, m);
            Z{i}(offset, 1) = 1;
            Z{i}(offset, 2) = pcpn(index, m);
            Z{i}(offset, 3) = zndx(index, m);
        end
    end
end
save('real_data/climate.mat', 'X', 'Z', 'Y', 'division_index');

real_data_split;