clear all
close all
clc

disp(' ')
filepath = '/Home/siv29/dsa030/Desktop/GEANT4/COSMIC_RAY_THUNDERSTORM-Geant4/build/timing_results.txt';
% yy = load(filepath);

fid = fopen(filepath);
tline = fgetl(fid);
read_vals = [];
while ischar(tline)
    read_vals(end+1,:) = str2num(tline(1:31));
    tline = fgetl(fid);
end

potentials = read_vals(:,1);
times = read_vals(:,6);

max_time_not_overflow = round(max(times(times<=20000-1)));



pot_list = unique(potentials);


for ii=1:length(pot_list)
    
    stats(ii) = sum(potentials==pot_list(ii));
    
    times_for_pot = times(potentials==pot_list(ii));
    
    cd_overflow = times_for_pot>=10000-5;
    nb_overflow(ii) = sum(cd_overflow);
    
    avg_run_time(ii) = mean(times_for_pot(~cd_overflow));
    
end

disp(num2str([pot_list(:) stats(:) avg_run_time(:) nb_overflow(:)]))

disp(' ')

floor(max(avg_run_time)*2)

disp(['Max simulation time (not overflow): ' num2str(max_time_not_overflow)])
disp(' ')