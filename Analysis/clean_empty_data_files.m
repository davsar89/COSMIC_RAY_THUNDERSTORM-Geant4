function clean_empty_data_files()

disp('cleaning empty datafiles...')

sett = load_settings();

IS_FRAM = sett.IS_FRAM;
base_path = sett.base_path;

POTENTIAL_LIST = sett.POTENTIAL_LIST;
RECORD_POS_LIST = sett.RECORD_POS_LIST;
EFIELD_CENTER_list = sett.EFIELD_CENTER_list;
EFIELD_SIZE_list = sett.EFIELD_SIZE_list;

for i_pot =1:length(POTENTIAL_LIST)
    for i_recPos = 1:length(RECORD_POS_LIST)
        for i_efield_c = 1:length(EFIELD_CENTER_list)
            for i_efield_s = 1:length(EFIELD_SIZE_list)
                
                POTENTIAL = POTENTIAL_LIST(i_pot);
                SELECTED_RECORD_POSITION = RECORD_POS_LIST(i_recPos);
                EFIELD_CENTER = EFIELD_CENTER_list(i_efield_c);
                EFIELD_SIZE = EFIELD_SIZE_list(i_efield_s);
                
                folder_simu_results_for_param_set = [ base_path num2str(POTENTIAL) 'MV/' num2str(SELECTED_RECORD_POSITION) '/' ...
                    num2str(EFIELD_CENTER) 'km_' num2str(EFIELD_SIZE) 'km/'];
                
                if ~exist(folder_simu_results_for_param_set, 'dir')
                    continue
                end
                
                data_files_for_run = dir([folder_simu_results_for_param_set '*/ALL_ener_mom_dists*.out']);
                
                for i_file = 1:length(data_files_for_run)
                    filename = [data_files_for_run(i_file).folder '/' data_files_for_run(i_file).name];
                    fid = fopen(filename);
                    if fgetl(fid) == -1
                        disp(['deleting : ' filename ])
                        fclose(fid);
                        delete(filename)
                        disp(['deleted: ' filename '   ;'])
                    else
                        fclose(fid);
                    end
                end
                
            end
        end
    end
end

disp('DONE.')

end
