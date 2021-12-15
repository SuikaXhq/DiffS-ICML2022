if is_demo==true
    fprintf('Generating data for demo case.\n');
    generate_data(50,5,100,10,10,0,4);
    return;
end

for case_number = 1:18

n = [100
    100
    100
    100
    100
    100
    100
    100
    100
    100
    100
    100
    100
    200
    300
    400
    500
    600];
M = [50
    100
    150
    200
    250
    300
    50
    50
    50
    50
    50
    50
    50
    50
    50
    50
    50
    50];
S = [5
    5
    5
    5
    5
    5
    3
    7
    9
    5
    5
    5
    5
    5
    5
    5
    5
    5];
p = [10
    10
    10
    10
    10
    10
    10
    10
    10
    20
    30
    40
    50
    10
    10
    10
    10
    10];
q = p;

n = n(case_number);
M = M(case_number);
S = S(case_number);
p = p(case_number);
q = q(case_number);


generate_data(M,S,n,p,q,case_number,10);

end
