clear all
close all
clc

sett=load_settings();

MIN_ENER_PHOT = 20; % keV
MAX_ENER_PHOT = 100000;
MIN_ENER_ELEC = 20;
MAX_ENER_ELEC = 100000;

idx_potential_0 = 4;

% ! sh COPY_FINAL_DATAFILES_FROM_FRAM.sh

if sett.IS_FRAM
    loaded = load('BIG_DATAFILE_all.mat');
else
    loaded = load([sett.base_path 'BIG_DATAFILE_all.mat']);
end

BIG_DATAFILE = loaded.BIG_DATAFILE;

POTENTIAL_LIST = sett.POTENTIAL_LIST;
RECORD_POS_LIST = sett.RECORD_POS_LIST;
EFIELD_CENTER_list = sett.EFIELD_CENTER_list;
EFIELD_SIZE_list = sett.EFIELD_SIZE_list;

ALT = 16;
RECORD_POS = 2;
EFIELD_SIZE = 1;

i_recPos = find_idx(RECORD_POS_LIST,RECORD_POS);
i_efield_c = find_idx(EFIELD_CENTER_list,ALT);
i_efield_s = find_idx(EFIELD_SIZE_list,EFIELD_SIZE);

nb_files=[];

%%

for i_pot =1:length(POTENTIAL_LIST)
    POTENTIAL = POTENTIAL_LIST(i_pot);
    
    if ~isempty(BIG_DATAFILE.photon{i_pot,i_recPos,i_efield_c,i_efield_s})
        nb_recorded(1,i_pot) = get_NB_RECORDED_in_energy_range(BIG_DATAFILE.photon{i_pot,i_recPos,i_efield_c,i_efield_s},MIN_ENER_PHOT,MAX_ENER_PHOT);
        nb_recorded(2,i_pot) = get_NB_RECORDED_in_energy_range(BIG_DATAFILE.electron{i_pot,i_recPos,i_efield_c,i_efield_s},MIN_ENER_ELEC,MAX_ENER_ELEC);
        nb_recorded(3,i_pot) = get_NB_RECORDED_in_energy_range(BIG_DATAFILE.positron{i_pot,i_recPos,i_efield_c,i_efield_s},MIN_ENER_ELEC,MAX_ENER_ELEC);
        
        nb_files(1,i_pot) = BIG_DATAFILE.photon{i_pot,i_recPos,i_efield_c,i_efield_s}.NB_FILES;
        nb_files(2,i_pot) = BIG_DATAFILE.electron{i_pot,i_recPos,i_efield_c,i_efield_s}.NB_FILES;
        nb_files(3,i_pot) = BIG_DATAFILE.positron{i_pot,i_recPos,i_efield_c,i_efield_s}.NB_FILES;
    else
        nb_recorded(1,i_pot) = 0;
        nb_recorded(2,i_pot) = 0;
        nb_recorded(3,i_pot) = 0;
    end
end

nb_recorded_before = nb_recorded;

for ii=1:size(nb_recorded,1)
    
    ref = nb_recorded(ii,6);
    if ref==0
        ref = interp1([-10 10],nb_recorded_before(ii,[idx_potential_0-1 idx_potential_0+1]),0);
    end
    nb_recorded(ii,:) = (nb_recorded(ii,:)-ref) ./ ref * 100.0 ;
    
end

disp(' ')
disp(num2str(nb_recorded))

disp(' ')
disp(num2str(nb_files))

% nb_recorded(nb_recorded>1000)=0;

%%
close all
plot(POTENTIAL_LIST,nb_recorded(1,:))
hold on
plot(POTENTIAL_LIST,nb_recorded(2,:))
plot(POTENTIAL_LIST,nb_recorded(3,:))
legend('photon','electron','positron')
grid on
xlabel('potential (MV)')
ylabel('Change w.r.t. 0 MV (%)')
rec_alt = ALT+RECORD_POS*EFIELD_SIZE/2.0;
title(['ALT = ' num2str(ALT) ' km; EFIELD SIZE = ' num2str(EFIELD_SIZE) ' km; RECORD POS = ' num2str(RECORD_POS) '   (=' num2str(rec_alt) ' km)'])
% set(gca,'yscale','log')

%%


function NB_RECORDED = get_NB_RECORDED_in_energy_range(data_struct,min_ener,max_ener)

grid = data_struct.ENERGY_GRID*1000.0;

eh = data_struct.ENERGY_HIST;

NB_RECORDED = sum(eh(grid>min_ener & grid<max_ener))*1e5;

end