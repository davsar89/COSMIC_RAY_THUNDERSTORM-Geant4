clear all
close all
clc

cd(fileparts(which(mfilename)));

% dbstop in ./clean_big_outlier_data_files at 76
% dbstop in ./clean_empty_data_files at 51

%%

sett=load_settings();

IS_FRAM = sett.IS_FRAM;

%% run all

% clean_empty_data_files();

% is_OK = false;
% while ~is_OK
%     try
%         clean_big_outlier_data_files();
%         is_OK=true;
%     catch ME
%         is_OK = false;
%     end
% end

try
    MAKE_BIG_DATAFILES;
catch ME
    try
        MAKE_BIG_DATAFILES;
    catch ME
        try
            MAKE_BIG_DATAFILES;
        catch ME
            try
                MAKE_BIG_DATAFILES;
            catch ME
                try
                    MAKE_BIG_DATAFILES;
                catch ME
                    try
                        MAKE_BIG_DATAFILES;
                    catch ME
                        try
                            MAKE_BIG_DATAFILES;
                        catch ME
                            try
                                MAKE_BIG_DATAFILES;
                            catch ME
                                try
                                    MAKE_BIG_DATAFILES;
                                catch ME
                                    MAKE_BIG_DATAFILES;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

clearvars -except IS_FRAM

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
                
                disp( [num2str([POTENTIAL SELECTED_RECORD_POSITION EFIELD_SIZE EFIELD_CENTER]) ' :: ' num2str(nb) ])
                
                total = total + nb;
                nb_sets = nb_sets+1;
            end
        end
    end
end

%%

cd(fileparts(which(mfilename)));

Generate_database_glow_simulations;

disp(['total: ' num2str(total)])

disp(['nb per param set: ' num2str(total/nb_sets)])

cd(fileparts(which(mfilename)));

%%
