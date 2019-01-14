function data = parse_output_file(filename)
data = [];
fid = fopen(filename);

if fgetl(fid) == -1
    disp(['deleting : ' filename ])
    delete(filename)
    disp('problem')
    error(['deleted: ' filename '   ; please re-run'])
else
    fclose(fid);
    fid = fopen(filename);
end

tline = fgetl(fid);

data.photon.ANGLE_ENERGY_HIST = zeros(13,256);
data.electron.ANGLE_ENERGY_HIST = zeros(13,256);
data.positron.ANGLE_ENERGY_HIST = zeros(13,256);

inc_rec_types = [0, 23, 23*2];

record_names={'photon','electron','positron'};

%   << settings->RANDOM_SEED
%   << settings->NB_EVENT
%   << settings->INITIAL_SAMPLE_TYPE
%   << settings->RECORD_TYPE
%   << settings->EFIELD_REGION_Y_CENTER
%   << settings->EFIELD_REGION_Y_FULL_LENGTH
%   << settings->POTENTIAL_VALUE
%   << settings->CURRENT_WEIGHT
%   << settings->RECORD_ALTITUDE
%   << settings->CR_GENERATION_ALT_MIN;

jj = 0;

for irt=1:length(record_names)
    i_fill_aeh=1;
    %disp(tline)
    for ii=1:21
        
        jj=jj+1;
        
        nums = str2num(tline);
        
        if jj==1+inc_rec_types(irt)
            data.(record_names{irt}).SEED = int32(nums(1));
            data.(record_names{irt}).SAMPLED_NB = nums(2);
            data.(record_names{irt}).INITIAL_SAMPLED_PDG = 777;
            data.(record_names{irt}).RECORD_PDG = nums(4);
            data.(record_names{irt}).EFIELD_ALT_CENTER = nums(5);
            data.(record_names{irt}).EFIELD_ALT_LENGTH = nums(6);
            data.(record_names{irt}).POTENTIAL = nums(7);
            data.(record_names{irt}).WEIGHT = nums(8);
            data.(record_names{irt}).RECORD_ALT = nums(9);
            data.(record_names{irt}).CR_INITIAL_ALT = nums(10);
        elseif jj==3+inc_rec_types(irt)
            data.(record_names{irt}).NB_RECORDED=nums(1);
        elseif jj==5+inc_rec_types(irt)
            data.(record_names{irt}).ZENITH_ANGLE_GRID=nums;
        elseif jj==7+inc_rec_types(irt)
            data.(record_names{irt}).ENERGY_GRID=nums;
        elseif jj>=9+inc_rec_types(irt)
            if ~isempty(nums)
                data.(record_names{irt}).ANGLE_ENERGY_HIST(i_fill_aeh,:) = nums;
                i_fill_aeh=i_fill_aeh+1;
            end
        end
        
        tline = fgetl(fid);
        
    end
end

% sett=load_settings();

% data.photon.ANGLE_ENERGY_HIST(data.photon.ANGLE_ENERGY_HIST>10^20)=sett.MIN_HIST_BCKGRND;
% data.electron.ANGLE_ENERGY_HIST(data.photon.ANGLE_ENERGY_HIST>10^20)=sett.MIN_HIST_BCKGRND;
% data.positron.ANGLE_ENERGY_HIST(data.photon.ANGLE_ENERGY_HIST>10^20)=sett.MIN_HIST_BCKGRND;
% 
% data.photon.ANGLE_ENERGY_HIST(data.photon.ANGLE_ENERGY_HIST<10^-20)=sett.MIN_HIST_BCKGRND;
% data.electron.ANGLE_ENERGY_HIST(data.photon.ANGLE_ENERGY_HIST<10^-20)=sett.MIN_HIST_BCKGRND;
% data.positron.ANGLE_ENERGY_HIST(data.photon.ANGLE_ENERGY_HIST<10^-20)=sett.MIN_HIST_BCKGRND;

fclose(fid);
end