clear all
close all
clc

set(0,'DefaultFigureWindowStyle','docked')
set(0, 'DefaultFigureRenderer', 'opengl');


cd(fileparts(which(mfilename)));

EVALUE_FLUX_INCREASE_ILDAS;
% EVALUE_FLUX_INCREASE_ALOFT1;
EVALUE_FLUX_INCREASE_EACK2000;
EVALUE_FLUX_INCREASE_EACK1996;
EVALUE_FLUX_INCREASE_ADELE2015;

cd(fileparts(which(mfilename)));

clear all
close all