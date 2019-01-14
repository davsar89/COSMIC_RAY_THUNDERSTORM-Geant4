clear all
close all
clc

yy = importdata('/Home/siv29/dsa030/Desktop/GEANT4/MOS_altitude/build/input/sampled_particles.txt');


type = yy(:,1);

ener = yy(:,2);

sum(type==1)
sum(type==2)
sum(type==3)

cos_angle = yy(:,3);
alt = yy(:,4);

figure(1)
hist(alt,120);

figure(2)
hist(-cos_angle,120);

figure(3)
make_spectrum(ener(type==1 & alt>11.9 & alt<12.1));
hold on
make_spectrum(ener(type==2 & alt>11.9 & alt<12.1));
make_spectrum(ener(type==3 & alt>11.9 & alt<12.1));

%% plot EXPACS data for check

plot_expcas(1)
plot_expcas(2)
plot_expcas(3)






%%


function [n_simu,xout] = make_spectrum(ener_list)
bins = logspace(-2,8,128);
[n_simu,xout] = histcounts(ener_list,bins);
n_simu = n_simu ./ diff(bins);
histogram('BinEdges',xout,'BinCounts',n_simu,'DisplayStyle','stairs','LineWidth',2);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
end


function plot_expcas(idx)

if (idx == 1)
    idx_part=9;
    factor = 200000;
elseif (idx == 2)
    idx_part=7;
    factor = 200000/3.9157e+03 * 15.9;
elseif (idx == 3)
    idx_part=8;
    factor = 200000/3.9157e+03;
end

filepath = '/Home/siv29/dsa030/Desktop/GEANT4/MOS_altitude/CR_data/12km_CR_spec_ildas_per_ener.txt';

fid = fopen(filepath);
C = textscan(fid,'%f %f %f %f %f %f %f %f %f','CommentStyle','#');
data = [C{:}];

eners = data(:,1);
phot_spec = data(:,idx_part);
energy_threshold=20;
index_kept = eners>energy_threshold/1000.;

eners = eners(index_kept);
phot_spec = phot_spec(index_kept);

hold on
loglog(eners,phot_spec/mean(phot_spec)*factor,'DisplayName','Initially sampled (EXPACS photon spectrum at 12 km)');

end