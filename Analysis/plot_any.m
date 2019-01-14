clear all
close all
clc

global POTENTIAL;
global SELECTED_RECORD_ALTITUDE;
global EFIELD_CENTER;
global EFIELD_SIZE;

POTENTIAL = 0;
SELECTED_RECORD_ALTITUDE = 10;
EFIELD_CENTER = 10;
EFIELD_SIZE = 1;

%bf = '/Home/siv29/dsa030/Desktop/GEANT4/COSMIC_RAY_THUNDERSTORM-Geant4/build/output/';
bf = '/Home/siv29/dsa030/Desktop/GEANT4/COSMIC_RAY_THUNDERSTORM-Geant4/build/output_FRAM/';
base_folder = [ bf num2str(POTENTIAL) 'MV/' num2str(SELECTED_RECORD_ALTITUDE) 'km/' ...
    num2str(EFIELD_CENTER) 'km_' num2str(EFIELD_SIZE) 'km/'];

it_is_balanced = check_if_balanced(base_folder);

if ~it_is_balanced
   disp('please re-balance data files')
%    return 
end

%% parsing all files
filelist = dir(fullfile(base_folder, '**/*.out'));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);

for ii=1:length(filelist)
    data{ii} = parse_output_file([filelist(ii).folder '/' filelist(ii).name]);
end

[N,BINS]=get_energy_histogram(data,22);
histogram('BinEdges',BINS,'BinCounts',N./diff(BINS),'displayStyle','stairs','LineWidth',2)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Energy (MeV)')
ylabel('Counts per bin (a.u.)')
hold on
[N,BINS]=get_energy_histogram(data,11);
histogram('BinEdges',BINS,'BinCounts',N./diff(BINS),'displayStyle','stairs','LineWidth',2)
[N,BINS]=get_energy_histogram(data,-11);
histogram('BinEdges',BINS,'BinCounts',N./diff(BINS),'displayStyle','stairs','LineWidth',2)
legend('photons','electrons','positrons')
grid on

[nb(1),nb_files_read1] = get_number_recorded(data,  22);
[nb(2),nb_files_read2] = get_number_recorded(data,  11);
[nb(3),nb_files_read3] = get_number_recorded(data, -11);

if length(unique([nb_files_read1 nb_files_read2 nb_files_read3]))~=1
    disp('error: data files were not re-balanced')
    return
end

scale=sum(nb)/nb_files_read1;
nb=nb./sum(nb)*100.0;


%% functions

function ener_grid = load_energy_grid()
ener_grid= [logspace(log10(0.003),log10(100.000),255) 1000000];
end

function [out_hist,bins] = get_energy_histogram(data, RECORD_PDG)
global SELECTED_RECORD_ALTITUDE;
bins = load_energy_grid();
out_hist=zeros(1,size(data{1}.ANGLE_ENERGY_HIST,2));
for ii=1:length(data)
    
    if data{ii}.RECORD_PDG == RECORD_PDG && data{ii}.RECORD_ALT == SELECTED_RECORD_ALTITUDE
        out_hist = out_hist + sum(data{ii}.ANGLE_ENERGY_HIST,1).*data{ii}.WEIGHT./double(data{ii}.SAMPLED_NB);
    end
    
end
out_hist=out_hist(1:end-1);
end

function [nb,nb_files_read] = get_number_recorded(data, RECORD_PDG)
global SELECTED_RECORD_ALTITUDE;
nb = 0;
nb_files_read = 0;
for ii=1:length(data)
    if data{ii}.RECORD_PDG == RECORD_PDG && data{ii}.RECORD_ALT == SELECTED_RECORD_ALTITUDE
        nb = nb + double(data{ii}.NB_RECORDED).*data{ii}.WEIGHT./double(data{ii}.SAMPLED_NB);
        nb_files_read = nb_files_read+1;
    end
    
end
end
