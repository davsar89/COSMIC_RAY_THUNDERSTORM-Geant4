function st=load_settings()

st={};

st.SMOOTHING_LEVEL = 0;

st.POTENTIAL_LIST = [-220 -210 -190 -175 -150 -120 -80 -50 -40 -30 -20 -10  0  10  20  30  40  50 80 120 150 175 190 210 220];
% st.POTENTIAL_LIST = [ -210 -190 -175 -150 -120 -80 -50 -40 -30 -20 -10  0  10  20  30  40  50 80 120 150 175 190 210];
st.RECORD_POS_LIST = [0  1 -1  2 -2 3 -3 4 -4  5 -5 6 -6  7 -7 8 -8];
st.EFIELD_CENTER_list = [6 8 10 12 14 15 16];
% st.EFIELD_CENTER_list = [6 8 10 12 14 16];
st.EFIELD_SIZE_list = [1 2 4];

st.POTENTIAL_LIMITS_FOR_CENTER_ALTITUDES = [310 250 200 160 125 90];

st.WANTED_RECORD_ALTS = [12, 14, 15, 20];

st.record_type_PDG_NB = [22, 11, -11];

st.record_names={'photon','electron','positron'};

st.MIN_HIST_BCKGRND = 10^-6;

[~, name] = system('hostname');

if contains(name,'iftrom039200')
    st.IS_FRAM = 0;
elseif contains(name,'fram.sigma2')
    st.IS_FRAM = 1;
else
    st.IS_FRAM = 1;
end

if st.IS_FRAM
    st.base_path = '/cluster/work/users/dsarria/SIMULATION_DATAFILES/COSMIC_THUNDER/';
else
    st.base_path = '/Data/ift/ift_romfys1/dsarria/SIMULATION_DATA/GLOW/';
end

st.i_0MV = find_idx(st.POTENTIAL_LIST,0);

st.POTENTIAL_LIST = sort(st.POTENTIAL_LIST);
st.EFIELD_CENTER_list = sort(st.EFIELD_CENTER_list);

st.responses_fodler = '/Home/siv29/dsa030/Desktop/article_cr_modeling/MATLAB_ANALYSIS/responses/';

st.figure_font_size = 14;

st.NB_SUPER_SAMPLING_PLOT=400;

st.particles_extra_weight=[1,1.0,1.0];

end