clear all
close all
clc

global POTENTIAL;
global SELECTED_RECORD_ALTITUDE;
global EFIELD_CENTER;
global EFIELD_SIZE;

POTENTIAL = 0;
SELECTED_RECORD_ALTITUDE = 16;
EFIELD_CENTER = 16;
EFIELD_SIZE = 1;

POTENTIAL_LIST = [0 5 10 15 20 25 30 35 40 45 50];
RECORD_ALTITUDE_LIST = [4 6 8 10 12 14 16];
EFIELD_CENTER_list = [4 6 8 10 12 14 16];
EFIELD_SIZE_list = [1];

%bf = '/Home/siv29/dsa030/Desktop/GEANT4/COSMIC_RAY_THUNDERSTORM-Geant4/build/output/';
bf = '/Home/siv29/dsa030/Desktop/GEANT4/COSMIC_RAY_THUNDERSTORM-Geant4/build/output_FRAM/';

%% Necessary to rebalance dataset because sometimes a given configuration buged and did not produce data
% e.g. there can be 5 files for initial photon but 4 files for initial proton; then one initial photon file should be removed

for i_pot =1:length(POTENTIAL_LIST)
    
    POTENTIAL = POTENTIAL_LIST(i_pot);
    
    for i_efield_s = 1:length(EFIELD_SIZE_list)
        for i_efield_c = 1:length(EFIELD_CENTER_list)
            
            EFIELD_SIZE = EFIELD_SIZE_list(i_efield_s);
            EFIELD_CENTER = EFIELD_CENTER_list(i_efield_c);
            
            for i_alt = 1:length(RECORD_ALTITUDE_LIST)
                
                SELECTED_RECORD_ALTITUDE = RECORD_ALTITUDE_LIST(i_alt);
                
                base_folder = [ bf num2str(POTENTIAL) 'MV/' num2str(SELECTED_RECORD_ALTITUDE) 'km/' ...
                    num2str(EFIELD_CENTER) 'km_' num2str(EFIELD_SIZE) 'km/'];
                
                if ~exist(base_folder, 'dir')
                    %disp('skipping since folder does not exist')
                    continue
                end
                
                disp(['POTENTIAL: ' num2str(POTENTIAL)])
                disp(['ALTITUDE: ' num2str(SELECTED_RECORD_ALTITUDE)])
                
                it_is_balanced=false;
                
                folders_types = {'initial_photon','initial_electron','initial_positron',...
                    'initial_muonN','initial_muonP','initial_neutron','initial_proton'};
                
                while ~it_is_balanced
                    folder_with_removed_stuff = [bf '/REMOVED'];
                    if ~exist(folder_with_removed_stuff, 'dir')
                        mkdir(folder_with_removed_stuff);
                    end
                    
                    %% parsing all files
                    filelist = dir(fullfile(base_folder, '**/*.out'));  %get list of files and folders in any subfolder
                    filelist = filelist(~[filelist.isdir]);
                    
                    to_make = [folder_with_removed_stuff '/empty' ];
                    if ~exist(to_make, 'dir')
                        mkdir(to_make)
                    end
                    
                    %% photons
                    the_folder = fullfile(base_folder, folders_types{1});
                    dir_list = dir(the_folder);
                    dir_list = dir_list([dir_list.isdir]);
                    dir_list = {dir_list([dir_list.isdir]).name};
                    dir_list = dir_list(~ismember(dir_list ,{'.','..'}));
                    nb_phot=0;
                    for ii=1:length(dir_list)
                        a=dir([the_folder '/' dir_list{ii} '/*.out']);
                        nb_out_files=size(a,1);
                        if nb_out_files==5
                            nb_phot=nb_phot+1;
                        else
                            movefile([the_folder '/' dir_list{ii}], [folder_with_removed_stuff '/empty/']);
                        end
                    end
                    
                    %% electrons
                    the_folder = fullfile(base_folder, folders_types{2});
                    dir_list = dir(the_folder);
                    dir_list = dir_list([dir_list.isdir]);
                    dir_list = {dir_list([dir_list.isdir]).name};
                    dir_list = dir_list(~ismember(dir_list ,{'.','..'}));
                    nb_elec=0;
                    for ii=1:length(dir_list)
                        a=dir([the_folder '/' dir_list{ii} '/*.out']);
                        nb_out_files=size(a,1);
                        if nb_out_files==5
                            nb_elec=nb_elec+1;
                        else
                            movefile([the_folder '/' dir_list{ii}], [folder_with_removed_stuff '/empty/']);
                        end
                    end
                    
                    %% positrons
                    the_folder = fullfile(base_folder, folders_types{3});
                    dir_list = dir(the_folder);
                    dir_list = dir_list([dir_list.isdir]);
                    dir_list = {dir_list([dir_list.isdir]).name};
                    dir_list = dir_list(~ismember(dir_list ,{'.','..'}));
                    nb_posi=0;
                    for ii=1:length(dir_list)
                        a=dir([the_folder '/' dir_list{ii} '/*.out']);
                        nb_out_files=size(a,1);
                        if nb_out_files==5
                            nb_posi=nb_posi+1;
                        else
                            movefile([the_folder '/' dir_list{ii}], [folder_with_removed_stuff '/empty/']);
                        end
                    end
                    
                    %% muon_n
                    the_folder = fullfile(base_folder, folders_types{4});
                    dir_list = dir(the_folder);
                    dir_list = dir_list([dir_list.isdir]);
                    dir_list = {dir_list([dir_list.isdir]).name};
                    dir_list = dir_list(~ismember(dir_list ,{'.','..'}));
                    nb_mun=0;
                    for ii=1:length(dir_list)
                        a=dir([the_folder '/' dir_list{ii} '/*.out']);
                        nb_out_files=size(a,1);
                        if nb_out_files==5
                            nb_mun=nb_mun+1;
                        else
                            movefile([the_folder '/' dir_list{ii}], [folder_with_removed_stuff '/empty/']);
                        end
                    end
                    
                    %% muon_p
                    the_folder = fullfile(base_folder, folders_types{5});
                    dir_list = dir(the_folder);
                    dir_list = dir_list([dir_list.isdir]);
                    dir_list = {dir_list([dir_list.isdir]).name};
                    dir_list = dir_list(~ismember(dir_list ,{'.','..'}));
                    nb_mup=0;
                    for ii=1:length(dir_list)
                        a=dir([the_folder '/' dir_list{ii} '/*.out']);
                        nb_out_files=size(a,1);
                        if nb_out_files==5
                            nb_mup=nb_mup+1;
                        else
                            movefile([the_folder '/' dir_list{ii}], [folder_with_removed_stuff '/empty/']);
                        end
                    end
                    
                    %% neutron
                    the_folder = fullfile(base_folder, folders_types{6});
                    dir_list = dir(the_folder);
                    dir_list = dir_list([dir_list.isdir]);
                    dir_list = {dir_list([dir_list.isdir]).name};
                    dir_list = dir_list(~ismember(dir_list ,{'.','..'}));
                    nb_neu=0;
                    for ii=1:length(dir_list)
                        a=dir([the_folder '/' dir_list{ii} '/*.out']);
                        nb_out_files=size(a,1);
                        if nb_out_files==5
                            nb_neu=nb_neu+1;
                        else
                            movefile([the_folder '/' dir_list{ii}], [folder_with_removed_stuff '/empty/']);
                        end
                    end
                    
                    %% proton
                    the_folder = fullfile(base_folder, folders_types{7});
                    dir_list = dir(the_folder);
                    dir_list = dir_list([dir_list.isdir]);
                    dir_list = {dir_list([dir_list.isdir]).name};
                    dir_list = dir_list(~ismember(dir_list ,{'.','..'}));
                    nb_pro=0;
                    for ii=1:length(dir_list)
                        a=dir([the_folder '/' dir_list{ii} '/*.out']);
                        nb_out_files=size(a,1);
                        if nb_out_files==5
                            nb_pro=nb_pro+1;
                        else
                            movefile([the_folder '/' dir_list{ii}], [folder_with_removed_stuff '/empty/']);
                        end
                    end
                    
                    %% rr-balance output files so that there is the same number for each particle type run
                    arr = [nb_phot nb_elec nb_posi nb_mun nb_mup nb_neu nb_pro];
                    
                    disp(num2str(arr))
                    
                    if length(unique(arr))==1
                        it_is_balanced = true;
                        disp('nothing to do')
                        continue
                    else
                        disp('necessary to re-balance data')
                    end
                    
                    [val,i_min]=min(arr);
                    
                    for ii=1:length(arr)
                        nb_to_remove(ii)  = arr(ii)-val;
                    end
                    
                    for i_type=1:7
                        
                        to_make = [folder_with_removed_stuff '/' folders_types{i_type}];
                        if ~exist(to_make, 'dir')
                            mkdir(to_make)
                        end
                        
                        the_folder = fullfile(base_folder, folders_types{i_type});
                        dir_list = dir(the_folder);
                        dir_list = dir_list([dir_list.isdir]);
                        dir_list = {dir_list([dir_list.isdir]).name};
                        dir_list = dir_list(~ismember(dir_list ,{'.','..'}));
                        
                        for dummy = 1:nb_to_remove(i_type)
                            idx_random_folder_to_remove = dummy;
                            
                            folder_if_already_moved = [folder_with_removed_stuff '/' folders_types{i_type} '/' dir_list{idx_random_folder_to_remove}];
                            if exist(folder_if_already_moved, 'dir')
                                system(['rm -rf ' folder_if_already_moved]);
                            end
                            movefile([the_folder '/' dir_list{idx_random_folder_to_remove}], ...
                                [folder_with_removed_stuff '/' folders_types{i_type} '/']);
                        end
                        
                    end
                    
                    it_is_balanced=check_if_balanced(base_folder);
                    if it_is_balanced
                        disp('Done. nothing more to do')
                        continue
                    end
                    
                end
            end
        end
    end
end