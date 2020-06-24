% Path and file specifications
pathname = '~/git_environment/Bayesian-matrixmodel/hunter-cevera/zinser/data/';
file = append(pathname, 'Zinser_SizeDist_calibrated-52-12.nc'); % Zinser data
count_file = append(pathname, 'Zinser_Figure2A.csv'); % Cell counts

num_days = 2;
if (num_days ~= 2) && (num_days ~= 1)
    error('Invalid number of days.')
end

% Extract all of the attributes of the Zinser data
time = double(ncread(file, 'time') + 60.0);
size_bounds = ncread(file, 'size_bounds');
w_obs = ncread(file, 'w_obs');
PAR = ncread(file, 'PAR');
m = ncread(file,'m');
delta_v_inv = ncread(file, 'delta_v_inv');
v_min = ncread(file, 'v_min');

% Read in cell counts
fig2 = readmatrix(count_file);
cell_counts = interp1(1:2:(num_days*24+1),fig2(:,2),time/60.0);
% cell_counts = fig2(:,2); % Cell culture A
% cell_counts = fig2(:,3); % Cell culture B
num_obs = num_days*12+1;
counts_as_matrix = zeros(m,num_obs);
for i = 1:num_obs
    counts_as_matrix(:,i) = cell_counts(i);
end

% Process data for input to the model
Edata = zeros(num_obs,2);
Edata(:,2) = PAR(1:num_obs);
% Edata(:,1) = time(1:num_obs)/60.0;
Edata(:,1) = int64(time/60.0);

volbins = size_bounds(1:m).';

N_dist = round(w_obs(1:num_obs,:).' .* counts_as_matrix);

save Zinser52_data Edata volbins N_dist;
