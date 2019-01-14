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

nb_run = 2

POTENTIAL_LIST = [-220, 220]

# if local run
if ("iftrom" in computer_name) or ("7370" in computer_name):
    POTENTIAL_LIST = [0]

REC_POS_list = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8]
Efield_alt_center_list = [6, 8, 10, 12, 14, 15, 16]
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
if ("iftrom" in computer_name) or ("7370" in computer_name):  # local (personal) computer

    nb_thread = 1  # number of threads (cpu) to run

    # Making an array where each element is the list of command for a given thread

    command_number = len(commands)

    print('Number of commands required ' + str(command_number))

    pool = Pool(nb_thread)  # to be always set to 1 for this MPI case
    for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands)):
        if returncode != 0:
            print("%d command failed: %d" % (i, returncode))

# COMPUTER CLUSTER RUN (uses MPI)
else:
    # MPI initializations and preliminaries
    comm = MPI.COMM_WORLD   # get MPI communicator object
    size = comm.Get_size()       # total number of processes
    rank = comm.Get_rank()       # rank of this process
    status = MPI.Status()   # get MPI status object

    # Define MPI message tags
    tags = enum('READY', 'DONE', 'EXIT', 'START')

    if rank == 0:
        # Master process executes code below
        tasks = commands
        task_index = 0
        num_workers = size - 1
        closed_workers = 0
        print("Master starting with %d workers" % num_workers)
        while closed_workers < num_workers:
            data = comm.recv(source=MPI.ANY_SOURCE,
                             tag=MPI.ANY_TAG, status=status)
            source = status.Get_source()
            tag = status.Get_tag()
            if tag == tags.READY:
                # Worker is ready, so send it a task
                if task_index < len(tasks):
                    comm.send(tasks[task_index], dest=source, tag=tags.START)
                    print("Sending task %d to worker %d" %
                          (task_index, source))
                    task_index += 1
                else:
                    comm.send(None, dest=source, tag=tags.EXIT)
            elif tag == tags.DONE:
                results = data
                print("Got data from worker %d" % source)
            elif tag == tags.EXIT:
                print("Worker %d exited." % source)
                closed_workers += 1

        print("Master finishing")
    else:
        # Worker processes execute code below
        name = MPI.Get_processor_name()
        print("I am a worker with rank %d on %s." % (rank, name))
        while True:
            comm.send(None, dest=0, tag=tags.READY)
            task = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
            tag = status.Get_tag()

            if tag == tags.START:
                # Do the work here
                task2 = [task]
                pool = Pool(1)  # to be always set to 1 for this MPI case
                for i, returncode in enumerate(pool.imap(partial(call, shell=True), task2)):
                    if returncode != 0:
                        print("%d command failed: %d" % (i, returncode))
                comm.send(returncode, dest=0, tag=tags.DONE)
            elif tag == tags.EXIT:
                break

    comm.send(None, dest=0, tag=tags.EXIT)
