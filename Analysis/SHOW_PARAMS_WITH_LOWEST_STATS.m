clear all
close all
clc

! sh COPY_FINAL_DATAFILES_FROM_FRAM.sh

sett=load_settings();

NB_TO_GET = 10;

loaded = load('BIG_DATAFILE_all.mat');
BIG_DATAFILE = loaded.BIG_DATAFILE;

POTENTIAL_LIST = sett.POTENTIAL_LIST;
RECORD_POS_LIST = sett.RECORD_POS_LIST;
EFIELD_CENTER_list = sett.EFIELD_CENTER_list;
EFIELD_SIZE_list = sett.EFIELD_SIZE_list;

NB_TOT = length(POTENTIAL_LIST)*length(EFIELD_SIZE_list)*length(EFIELD_CENTER_list)*length(RECORD_POS_LIST);

nb_err = 0;

for i_pot = 1:length(POTENTIAL_LIST)
    for i_efield_s = 1:length(EFIELD_SIZE_list)
        for i_efield_c = 1:length(EFIELD_CENTER_list)
            for i_recPos = 1:length(RECORD_POS_LIST)
                
                if ~isempty(BIG_DATAFILE.photon{i_pot,i_recPos,i_efield_c,i_efield_s})
                    stats(i_pot,i_recPos,i_efield_c,i_efield_s) = BIG_DATAFILE.photon{i_pot,i_recPos,i_efield_c,i_efield_s}.NB_FILES;
                    
                    nb_f_phot = BIG_DATAFILE.photon{i_pot,i_recPos,i_efield_c,i_efield_s}.NB_FILES;
                    nb_f_elec = BIG_DATAFILE.electron{i_pot,i_recPos,i_efield_c,i_efield_s}.NB_FILES;
                    nb_f_posi = BIG_DATAFILE.positron{i_pot,i_recPos,i_efield_c,i_efield_s}.NB_FILES;
                    
                    if (nb_f_phot~=nb_f_elec)
                        disp('nb_f_phot should be the same as nb_f_elec')
                        nb_err = nb_err+1;
                    elseif (nb_f_elec~=nb_f_posi)
                        disp('nb_f_elec should be the same as nb_f_posi')
                        nb_err = nb_err+1;
                    end
                    
                else
                    stats(i_pot,i_recPos,i_efield_c,i_efield_s) = 0;
                end
            end
        end
    end
end

disp(['Number with 0: ' num2str(sum(stats(:)==0))  ' / ' num2str(NB_TOT)])

disp(['nb of inconsistent photon/electron/positron number of files processed (probably because simulations are still running):' num2str(nb_err) ' / ' num2str(NB_TOT)])

min_stats = 10.e6;

i_found = 1;

while i_found<=NB_TO_GET
    
    for i_pot =1:length(POTENTIAL_LIST)
        for i_recPos = 1:length(RECORD_POS_LIST)
            for i_efield_c = 1:length(EFIELD_CENTER_list)
                for i_efield_s = 1:length(EFIELD_SIZE_list)
                    
                    if stats(i_pot,i_recPos,i_efield_c,i_efield_s)<min_stats
                        
                        lowest_idx_list{i_found} = [i_pot,i_recPos,i_efield_c,i_efield_s];
                        
                        min_stats = stats(i_pot,i_recPos,i_efield_c,i_efield_s);
                        
                        stats_list(i_found) = min_stats;
                        
                        stats(i_pot,i_recPos,i_efield_c,i_efield_s) = 10.e6;
                        
                    end
                    
                end
            end
        end
    end
    
    i_found = i_found+1;
    min_stats = 10.e6;
    
end


for ii=1:length(stats_list)
    
    pot = POTENTIAL_LIST(lowest_idx_list{ii}(1));
    recPos = RECORD_POS_LIST(lowest_idx_list{ii}(2));
    efield_c = EFIELD_CENTER_list(lowest_idx_list{ii}(3));
    efield_s = EFIELD_SIZE_list(lowest_idx_list{ii}(4));
    
    disp([ num2str([pot recPos efield_c efield_s]) ' :: ' num2str(stats_list(ii))])
    
end
