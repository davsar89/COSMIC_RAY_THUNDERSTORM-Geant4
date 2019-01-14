clear all
close all
clc

%%

response=[];
response.name='ILDAS';
response.energy_grid =             [30  40  50  60  70  80  90  100  110  120 200 300 400 500 600 700 800 900 1000  10000]; % keV
response.effective_area.photon =   [10  20  30  40  50  50  50  50   50   50  50  60  60  70  70  70  80  80  90      110]; % cm2
response.effective_area.electron = [ 0  00  00  00  00  00  00  10   10   20  20  20  20  20  20  20  30  30  30       50]; % cm2
response.effective_area.positron = [ 0  00  00  00  00  00  00  10   10   20  20  20  20  20  20  20  30  30  30       50]; % cm2
response.effective_area.unit='cm2';
response.energy_unit='keV';

responses_fodler = '/Home/siv29/dsa030/Desktop/article_cr_modeling/MATLAB_ANALYSIS/responses/';

save([responses_fodler 'EACK1996.mat'],'response');