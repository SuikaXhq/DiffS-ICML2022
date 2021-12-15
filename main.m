is_demo = true;
% is_demo = false;
data_generate;
calc_unit_GLS;
main_full('dishes', is_demo);
main_full('kmeans', is_demo);
main_full('kmeans-open', is_demo);
main_full('cd', is_demo);

