import platform
import subprocess as sp
import time
import re
from os import system
from os import walk
import numpy as np
from os import listdir
from mpi4py import MPI
from subprocess import call
from functools import partial
from multiprocessing.dummy import Pool
import numpy as np
import time
import random
import sys
from random import randint
from time import sleep

computer_name = platform.node()

# functions definitions for MPI usage (optional)


def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

################################################################################
# Defining commands to run

#sleep(randint(1, 5))

nb_run = 2

#POTENTIAL_LIST= [0]
#POTENTIAL_LIST= [0, 10, 20, 30, 40, 50]
#POTENTIAL_LIST= [50, 40, 30, 20, 10, 0]
#POTENTIAL_LIST = [-120, -80, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 80, 120, -150, 150, -175, 175, -190, 190]
#POTENTIAL_LIST = [-210, 210, -220, 220]
POTENTIAL_LIST = [-220, -210, -190, -175, -150, -120, -80, -50, -40, -30, -20, -10, 
                    0,
                    10, 20, 30, 40, 50, 80, 120, 150, 175, 190, 210, 220]
#POTENTIAL_LIST = [-150, 150, -175, 175, -190, 190]
#POTENTIAL_LIST = [-50, -40, -30, -20, -10]

# if local run
#if ("iftrom" in computer_name) or ("7370" in computer_name):
#    POTENTIAL_LIST = [0]

#REC_POS_list=           [0]
REC_POS_list = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8]
#REC_POS_list = [2]
#Efield_alt_center=      [16]
Efield_alt_center_list = [4, 6, 8, 10, 11, 12, 13, 14, 15, 16]
#Efield_alt_center_list = [10]
#Efield_alt_center=      [12,  14]
#Efield_full_size=       [1, 2]
#Efield_full_size=       [1, 2, 4]
Efield_full_size_list = [2]

# defining the commands to be run in parallel
commands = []
executable = 'timeout 20000s ./mos_test'

for _ in range(nb_run):
    for ii, _ in enumerate(Efield_alt_center_list):
        for jj, _ in enumerate(REC_POS_list):
            for kk, _ in enumerate(Efield_full_size_list):
                for POT in POTENTIAL_LIST:

                    rec_alt = Efield_alt_center_list[ii] + REC_POS_list[jj]*Efield_full_size_list[kk]/2.0

                    if rec_alt not in [12, 14, 15, 20]:
                        continue

                    if abs(POT) == 0:
                        INITIAL_NUMBER_TO_SAMPLE = 1000000
                    elif abs(POT) == 10:
                        INITIAL_NUMBER_TO_SAMPLE = 50000
                    elif abs(POT) == 20:
                        INITIAL_NUMBER_TO_SAMPLE = 50000
                    elif abs(POT) == 30:
                        INITIAL_NUMBER_TO_SAMPLE = 50000
                    elif abs(POT) == 40:
                        INITIAL_NUMBER_TO_SAMPLE = 50000
                    elif abs(POT) == 50:
                        INITIAL_NUMBER_TO_SAMPLE = 50000
                    elif abs(POT) == 80:
                        INITIAL_NUMBER_TO_SAMPLE = 20000
                    elif abs(POT) == 120:
                        INITIAL_NUMBER_TO_SAMPLE = 20000
                    elif abs(POT) > 120:
                        INITIAL_NUMBER_TO_SAMPLE = 10000

                    commands.append(
                        executable
                        + ' ' + str(POT)
                        + ' ' + str(INITIAL_NUMBER_TO_SAMPLE)
                        + ' ' + str(REC_POS_list[jj])
                        + ' ' + str(Efield_alt_center_list[ii])
                        + ' ' + str(Efield_full_size_list[kk])
                    )

#####################################

random.shuffle(commands)

# LOCAL RUN (uses python multiprocessing library)
nb_thread = 4  # number of threads (cpu) to run

# Making an array where each element is the list of command for a given thread

command_number = len(commands)

print('Number of commands required ' + str(command_number))

pool = Pool(nb_thread)  # to be always set to 1 for this MPI case
for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands)):
    if returncode != 0:
        print("%d command failed: %d" % (i, returncode))