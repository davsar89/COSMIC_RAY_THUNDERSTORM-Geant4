import platform
import subprocess as sp
import time
import re
import os
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
from datetime import datetime

computer_name = platform.node()

# functions definitions for MPI usage (optional)

def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

################################################################################
# Defining commands to run

nb_run = 10

POTENTIAL_LIST = [-220, -210, -190, -175, -150, -120, -80, -50, -40, -30, -20, -10, 
                    0,
                    10, 20, 30, 40, 50, 80, 120, 150, 175, 190, 210, 220]

REC_POS_list = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8]
Efield_alt_center_list = [6, 8, 10, 12, 14, 15, 16]
Efield_full_size_list = [2]

# defining the commands to be run in parallel
timelimit = 20000 # sec
commands = []
executable = f"timeout {timelimit}s ./mos_test"

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
                        INITIAL_NUMBER_TO_SAMPLE = 500000
                    elif abs(POT) == 20:
                        INITIAL_NUMBER_TO_SAMPLE = 100000
                    elif abs(POT) == 30:
                        INITIAL_NUMBER_TO_SAMPLE = 100000
                    elif abs(POT) == 40:
                        INITIAL_NUMBER_TO_SAMPLE = 500000
                    elif abs(POT) == 50:
                        INITIAL_NUMBER_TO_SAMPLE = 500000
                    elif abs(POT) == 80:
                        INITIAL_NUMBER_TO_SAMPLE = 200000
                    elif abs(POT) == 120:
                        INITIAL_NUMBER_TO_SAMPLE = 200000
                    elif abs(POT) > 120:
                        INITIAL_NUMBER_TO_SAMPLE = 100000

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

if __name__ == "__main__":

    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    #with open("timing_results.txt", "w+") as file_object:
    #    file_object.write('')

    for comm in commands:
        t0 = time.perf_counter()
        os.system(comm)
        t1 = time.perf_counter()
        ellapsed = t1-t0

        with open("timing_results.txt", "a") as file_object:
            now = datetime.now() # current date and time
            date_time_str = now.strftime("%m/%d/%Y, %H:%M:%S")
            file_object.write(comm.replace(f"{executable} ", '') + ' ' + str(ellapsed) + ' :: ' + date_time_str + '\n')
