clear all
close all
clc

disp(' ')

%% display stats status

sett=load_settings();

IS_FRAM = sett.IS_FRAM;

if ~IS_FRAM
    load([sett.base_path 'BIG_DATAFILE_all.mat'])
else
    load(['BIG_DATAFILE_all.mat'])
end

bdf = BIG_DATAFILE.photon;

sett=load_settings();
POTENTIAL_LIST = sett.POTENTIAL_LIST;
RECORD_POS_LIST = sett.RECORD_POS_LIST;
EFIELD_CENTER_list = sett.EFIELD_CENTER_list;
EFIELD_SIZE_list = sett.EFIELD_SIZE_list;

total = 0;

nb_sets = 0;

for i_pot =1:length(POTENTIAL_LIST)
    for i_efield_s = 1:length(EFIELD_SIZE_list)
        for i_efield_c = 1:length(EFIELD_CENTER_list)
            for i_alt = 1:length(RECORD_POS_LIST)
                
                POTENTIAL = POTENTIAL_LIST(i_pot);
                SELECTED_RECORD_POSITION = RECORD_POS_LIST(i_alt);
                EFIELD_SIZE = EFIELD_SIZE_list(i_efield_s);
                EFIELD_CENTER = EFIELD_CENTER_list(i_efield_c);
                
                dat = bdf{i_pot,i_alt,i_efield_c,i_efield_s};
                
                if isempty(dat)
                    nb=0;
                else
                    nb=dat.NB_FILES;
                end
                
                rec_alt = EFIELD_CENTER + SELECTED_RECORD_POSITION*EFIELD_SIZE/2.0;
                
                if ismember(rec_alt,sett.WANTED_RECORD_ALTS) && POTENTIAL<sett.POTENTIAL_LIMITS_FOR_CENTER_ALTITUDES(i_efield_c)
                    disp( [num2str([POTENTIAL SELECTED_RECORD_POSITION EFIELD_SIZE EFIELD_CENTER]) ' :: ' num2str(nb) ])
                end
                
                total = total + nb;
                nb_sets = nb_sets+1;
            end
        end
    end
end

%%

disp(['total: ' num2str(total)])

disp(['nb per param set: ' num2str(total/nb_sets)])