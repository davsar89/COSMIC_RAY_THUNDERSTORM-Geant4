function it_is_balanced=check_if_balanced(base_folder)

%% photons
the_folder = fullfile(base_folder, 'initial_photon');
dir_list = dir(the_folder);
dir_list = dir_list([dir_list.isdir]);
dir_list = {dir_list([dir_list.isdir]).name};
dir_list = dir_list(~ismember(dir_list ,{'.','..'}));
nb_phot=0;
for ii=1:length(dir_list)
    a=dir([the_folder '/' dir_list{ii} '/*.out']);
    nb_out_files=size(a,1);
    if nb_out_files==3
        nb_phot=nb_phot+1;
    end
end

%% electrons
the_folder = fullfile(base_folder, 'initial_electron');
dir_list = dir(the_folder);
dir_list = dir_list([dir_list.isdir]);
dir_list = {dir_list([dir_list.isdir]).name};
dir_list = dir_list(~ismember(dir_list ,{'.','..'}));
nb_elec=0;
for ii=1:length(dir_list)
    a=dir([the_folder '/' dir_list{ii} '/*.out']);
    nb_out_files=size(a,1);
    if nb_out_files==3
        nb_elec=nb_elec+1;
    end
end

%% positrons
the_folder = fullfile(base_folder, 'initial_positron');
dir_list = dir(the_folder);
dir_list = dir_list([dir_list.isdir]);
dir_list = {dir_list([dir_list.isdir]).name};
dir_list = dir_list(~ismember(dir_list ,{'.','..'}));
nb_posi=0;
for ii=1:length(dir_list)
    a=dir([the_folder '/' dir_list{ii} '/*.out']);
    nb_out_files=size(a,1);
    if nb_out_files==3
        nb_posi=nb_posi+1;
    end
end

%% muon_n
the_folder = fullfile(base_folder, 'initial_muonN');
dir_list = dir(the_folder);
dir_list = dir_list([dir_list.isdir]);
dir_list = {dir_list([dir_list.isdir]).name};
dir_list = dir_list(~ismember(dir_list ,{'.','..'}));
nb_mun=0;
for ii=1:length(dir_list)
    a=dir([the_folder '/' dir_list{ii} '/*.out']);
    nb_out_files=size(a,1);
    if nb_out_files==3
        nb_mun=nb_mun+1;
    end
end

%% muon_p
the_folder = fullfile(base_folder, 'initial_muonP');
dir_list = dir(the_folder);
dir_list = dir_list([dir_list.isdir]);
dir_list = {dir_list([dir_list.isdir]).name};
dir_list = dir_list(~ismember(dir_list ,{'.','..'}));
nb_mup=0;
for ii=1:length(dir_list)
    a=dir([the_folder '/' dir_list{ii} '/*.out']);
    nb_out_files=size(a,1);
    if nb_out_files==3
        nb_mup=nb_mup+1;
    end
end

%% neutron
the_folder = fullfile(base_folder, 'initial_neutron');
dir_list = dir(the_folder);
dir_list = dir_list([dir_list.isdir]);
dir_list = {dir_list([dir_list.isdir]).name};
dir_list = dir_list(~ismember(dir_list ,{'.','..'}));
nb_neu=0;
for ii=1:length(dir_list)
    a=dir([the_folder '/' dir_list{ii} '/*.out']);
    nb_out_files=size(a,1);
    if nb_out_files==3
        nb_neu=nb_neu+1;
    end
end

%% proton
the_folder = fullfile(base_folder, 'initial_proton');
dir_list = dir(the_folder);
dir_list = dir_list([dir_list.isdir]);
dir_list = {dir_list([dir_list.isdir]).name};
dir_list = dir_list(~ismember(dir_list ,{'.','..'}));
nb_pro=0;
for ii=1:length(dir_list)
    a=dir([the_folder '/' dir_list{ii} '/*.out']);
    nb_out_files=size(a,1);
    if nb_out_files==3
        nb_pro=nb_pro+1;
    end
end

%%

arr = [nb_phot nb_elec nb_posi nb_mun nb_mup nb_neu nb_pro];

disp(num2str(arr))
pause(0.5)

if sum(arr)==0
    it_is_balanced=false;
    disp('WARNING: there are no files')
end

if length(unique(arr))==1
    disp('it is balanced')
    it_is_balanced=true;
else
    it_is_balanced=false;
    disp('it is NOT balanced')
end


end