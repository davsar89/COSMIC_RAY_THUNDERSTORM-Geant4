import subprocess as sp
import time
import random
import re
from os import walk
import numpy as np
from os import listdir
from os.path import isfile, join, getsize


######

re1='(Particle)'	# Word 1
re2='(_)'	# Any Single Character 1
re3='(List)'	# Word 2
re4='(_)'	# Any Single Character 2
re5='(\\d+)'	# Integer Number 1
re6='(_)'	# Any Single Character 3
re7='(120)'	# Integer Number 2
re8='(km)'	# Word 3
re9='(_)'	# Any Single Character 4
re10='(-?)'	# Any Single Character 5
re11='(\\d+)'	# Integer Number 3
re12='(_)'	# Any Single Character 6

rg_potential_read = re.compile(re1+re2+re3+re4+re5+re6+re7+re8+re9+re10+re11+re12,re.IGNORECASE|re.DOTALL)

##########################################################
##########################################################

def get_batch_script_content(python_file,buid_folder,nb_cpu,potential):
  potential=int(potential)
  potential_str=str(potential)
  if (potential>=0):
      potential_str="+"+potential_str

  batch_script="""#!/bin/bash
## Project:
##SBATCH --account=NN9526K --partition=bigmem
#SBATCH --account=NN9526K --qos=preproc
## Job name:
#SBATCH --job-name=_{3}_
## Wall time limit:
#SBATCH --time=0-2:0:0
## Number of nodes and task per node
#SBATCH --nodes=1 --ntasks-per-node=1  --cpus-per-task={2}
#SBATCH --mail-user=david.sarria@uib.no
##SBATCH --mem-per-cpu=2G

## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## Software modules
module restore system   # Restore loaded modules to the default

module load foss/2017a
module load CMake/3.9.1
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.1-GCC-6.4.0-2.28
module load Python/3.6.4-foss-2018a

module list             # List loaded modules, for easier debugging

## Run the application
echo "Running code..."
cd {1}
srun python {0}
echo "Done."
""".format(python_file,buid_folder,nb_cpu,potential_str)
  
  return batch_script

##########################################################
##########################################################

def get_python_script_content(potential,nb_shoot,nb_cpu):
    python_script = """#!/usr/bin/env python2.7
import platform
from subprocess import call
from functools import partial
from multiprocessing.dummy import Pool
import numpy as np
import random
import time
import sys
from datetime import datetime
import subprocess as sp

computer_name = platform.node()

#### functions definitions

############################

random.seed(datetime.now())
seedd = random.randint(1,1000000000)

################################################################################
### Defining commands to run

#efield_altitudes =  [9.,10.,11.,12.,13.,14.]# kilometers, record altitude

efield_altitudes = [21.0]# kilometers, record altitude

scales = [1.5436, 1.9256, 2.4605, 3.2461, 4.4152, 6.0826, 8.3489, 11.4043]

efield_sizes = [2.0] # does not matter if efield_altitudes > 20

#potential_list = [40., 60., 80. ,100., 120., 140., 160.]
record_altitudes = [20.0]

# -1 and 1 are used as a proxy of zero field since zero field leads to a not fully understood behaviour.
potential_list_0 = np.array([-10.,10.,-30.,30.,-60.,60.,-80.,80.,-100.,100.,-150.,150.,-200.,200.,-300.,300.,-400.,400.]) # will be changed
#potential_list_0 = potential_list_0[potential_list_0<149.]

tilt_list = np.array([0.]) # does not matter if efield_altitudes > 20
#tilt_list = np.array([0.,25.,45.])

if ("iftrom" in computer_name) or ("7370" in computer_name) or ("sarria-pc" in computer_name):
    nb_run = 1
else :
    nb_run = 1000
#potential_list = [200.]

nb_record_to_shoot_per_run = {1}

# defining the commands to be run in parallel

commands=[]
excecutable = './mos_test'

#seedd += 1

for _ in range(nb_run):
    for size in efield_sizes:
        for alti_e in efield_altitudes:
        
            potential_list_0 = np.array([{0}])
            
            for pot in potential_list_0:
                for tilt in tilt_list:
                    commands.append(excecutable + ' ' + str(seedd) + ' ' + str(nb_record_to_shoot_per_run) 
                                    + ' ' + str(alti_e) + ' ' + str(size) + ' ' + str(pot) + ' ' + str(tilt)
                                    + ' ' + str(record_altitudes[0]))
                    seedd+=1

################################################################################
print(len(commands))
print(commands[0])
#############

nb_thread = int({2}) # number of threads (cpu) to run

commands2 = commands[0:int(nb_thread)-1] # one command per thread

command_number = len(commands2)

#print('Number of commands required '+ str(command_number))

pool = Pool(nb_thread) #
for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands2)):
    if returncode != 0:
        print("%d command failed: %d" % (i, returncode))

""".format(potential,nb_shoot,nb_cpu)
  
    return python_script
  
##############################################################################
##############################################################################

def get_potential(text_file):
  re1='(_)'  # Any Single Character 1
  re2='(\\d+)'  # Integer Number 1
  re3='(_)'  # Any Single Character 2
  re4='(\\d+)'  # Integer Number 2
  re5='(km_)'  # Variable Name 1
  re6='([-+]\\d+)'  # Integer Number 1
  re7='(_)'  # Any Single Character 3
  re8='(\\d+)'  # Integer Number 3
  re9='.*?'  # Non-greedy match on filler
  re10='(\\d+)'  # Integer Number 4
  re11='(\\.)'  # Any Single Character 4
  re12='(out)'  # Word 1

  rg = re.compile(re1+re2+re3+re4+re5+re6+re7+re8+re9+re10+re11+re12,re.IGNORECASE|re.DOTALL)
  mg_search = rg.search(text_file)
  #print(mg_search)
  
  try:
      int3 = mg_search.group(6)
  except:
    try:
      re6='(\\d+)'  # Integer Number 1
      rg = re.compile(re1+re2+re3+re4+re5+re6+re7+re8+re9+re10+re11+re12,re.IGNORECASE|re.DOTALL)
      mg_search = rg.search(text_file)
      int3 = mg_search.group(6)
    except:
      int3 = 123456789
      
  return int3


def make_running_id_and_potential_list():
    ## getting job ID list
    bashCommand = "squeue -u dsarria"
    process = sp.Popen(bashCommand.split(), stdout=sp.PIPE)
    output, error = process.communicate()
    output = output.splitlines()
    output = output[1:]
    #print(len(output))
    #print(output)
    
    re1='(\\d+)'	# Integer Number 1
    re2='(\\s+)'	# White Space 1
    re3='((?:[a-z][a-z]+))'	# Word 1
    re4='(\\s+)'	# White Space 2
    re5='(_)'	# Any Single Character 1
    re6='([-+]\\d+)'	# Integer Number 1
    re7='(_)'	# Any Single Character 2
    re8='(\\s+)'	# White Space 3
    re9='(dsarria)'	# Word 2
    re10='(\\s+)'	# White Space 4
    re11='((?:[a-z][a-z0-9_]*))'	# Variable Name 1
    re12='(\\s+)'	# White Space 5

    rg = re.compile(re1+re2+re3+re4+re5+re6+re7+re8+re9+re10+re11+re12,re.IGNORECASE|re.DOTALL)
    
    running_id_list = []
    running_pot_list = []
    
    for element in output:
      m = rg.search(element) 
      if m:
        slurm_process_id = m.group(1)
        potential = int(m.group(6))
        running_char = str(m.group(11))

        #print(running_char)

        if (running_char=="R") or (running_char=="PD"):
          running_id_list.append(int(slurm_process_id))
          running_pot_list.append(int(potential))
 
    return running_id_list, running_pot_list

##############################################################################
##############################################################################
  
def check_if_already_running(potential,ID_list,potential_list,max_nb_proc_per_pot):
    if ID_list==[]:
        return False
      
    ID_list, potential_list = make_running_id_and_potential_list()

    #print(potential_list)
    if int(potential) in potential_list:
        if (potential_list.count(int(potential)) >= max_nb_proc_per_pot):
            return True
    else:
        return False

##############################################################################
##############################################################################
      
def find_potential_from_filename(filename):

    m = rg_potential_read.search(filename)
    
    if m:
        c5   = m.group(10)
        int3 = int(m.group(11))
        
    if (c5=="-") :
        sign = -1
    else:
        sign = 1
       
    potential = sign*int3
     
    return potential

##############################################################################
##############################################################################

def check_statistics(potential_list):
  print("starting check_statistics...")
  mypath = "/cluster/work/users/dsarria/G4_projects/BELOW_new_method/build/output/"

  # getting list of all output files

  files = [f for f in listdir(mypath) if isfile(join(mypath, f))]

  print("output file list retrieved.")
  files = list(set(files)) # to remove duplicates, just in case

  #print(files)

  # checking the particle number recorded for each potential
  statistis_dict = {}

  # initilisation with zeros
  for pot in potential_list:
    key = str(int(float(pot)))
    statistis_dict[key] = 0
    
  kB_per_particle = 0.2338 # size (kiloBytes) for each line of record file

  # counting (may introduce new keys)
  total = len(files)
  done_idx = 0
  
  for file_name in files:
    size = float(getsize(mypath+file_name))/1000.0 # bytes to kbytes
    
    with open(mypath+file_name) as my_file:
      
        first_char = my_file.read(1) #get the first character to check if file is empty or not
        if first_char: # file in not empty
  
          #first_line = my_file.readline() # loading the first line

          #loaded = np.fromstring(first_line, sep=' ')
        
          #potential = loaded[16]
          potential = find_potential_from_filename(file_name)
          #print(potential)
          
          #num_lines = 0
          #for line in my_file.xreadlines( ):
              #num_lines += 1
          num_lines = int(size / kB_per_particle)
          
          #print(num_lines)

          key = str(int(float(potential)))
      
          if key in statistis_dict:
              statistis_dict[key] = statistis_dict[key] + num_lines
          else:
              statistis_dict[key] = num_lines

    done_idx += 1
    #print("Done {} / {}".format(done_idx,total))

  for key in statistis_dict:
      print("{} : {}".format(key,statistis_dict[key]))

  return statistis_dict


##############################################################################
##############################################################################

def check_statistics_v2(potential_list):
  print("starting check_statistics...")
  mypath = "/cluster/work/users/dsarria/G4_projects/BELOW_new_method/build/output/"
  
  statistis_dict = {}

  # initilisation with zeros
  for pot in potential_list:
      key = str(int(float(pot)))
      statistis_dict[key] = 0


  kB_per_particle = 0.2338*1000.0 # average size (Bytes) for each line of record file

  bashCommand = "ls -l " + mypath
  process = sp.Popen(bashCommand.split(), stdout=sp.PIPE)
  output2, error = process.communicate()
  output2 = output2.replace('\n',' ')
  
  #print(output2)

  if (output2 == "total 0"):
      return statistis_dict
  
  re1='(dsarria)'	# Word 1
  re2='(\\s+)'	# White Space 1
  re3='(\\d+)'	# Integer Number 1
  re4='(\\s+)'	# White Space 2
  re5='((?:[a-z][a-z0-9_]*))'	# Variable Name 1
  re6='(\\s+)'	# White Space 3
  re7='(\\d+)'	# Integer Number 2
  re8='(\\s+)'	# White Space 4
  re9='((?:(?:[0-1][0-9])|(?:[2][0-3])|(?:[0-9])):(?:[0-5][0-9])(?::[0-5][0-9])?(?:\\s?(?:am|AM|pm|PM))?)'
  re10='(\\s+)'	# White Space 5
  re11='(Particle)'	# Word 2
  re12='(_)'	# Any Single Character 1
  re13='(List)'	# Word 3
  re14='(_)'	# Any Single Character 2
  re15='(\\d+)'	# Integer Number 3
  re16='(_)'	# Any Single Character 3
  re17='(\\d+)'	# Integer Number 4
  re18='(km)'	# Word 4
  re19='(_)'	# Any Single Character 4
  re20='(-?)'	# Any Single Character 5
  re21='(\\d+)'	# Integer Number 5
  re22='(_)'	# Any Single Character 6

  rg = re.compile(re1+re2+re3+re4+re5+re6+re7+re8+re9+re10+re11+re12+re13+re14+re15+re16+re17+re18+re19+re20+re21+re22,re.IGNORECASE|re.DOTALL)

  found_list = rg.findall(output2)
  
  #print(found_list)

  #potential_list = []
  #size_list = []

  if found_list:
      for element in found_list:
        file_size = int(element[2])
        potential = int(element[20])
    
        possible_sign = element[19]
    
        if(possible_sign=="-"):
            sign=-1
        else:
            sign=1
    
        potential = potential*sign
    
        #potential_list.append(potential)
        #size_list.append(file_size)
        
        key = str(int(float(potential)))
        
        num_lines = int(float(file_size)/kB_per_particle)
      
        if key in statistis_dict:
            statistis_dict[key] = statistis_dict[key] + num_lines
        else:
            statistis_dict[key] = num_lines

  return statistis_dict


##############################################################################
##############################################################################

def launch_job(potentia, nb_shoo, nb_cpuu, a_random_num, iidx, process_ID_lst, potential_lst):
    max_nb_proc_per_pot = 3

    if (check_if_already_running(potentia,process_ID_lst,potential_lst,max_nb_proc_per_pot)):
      print("Not starting potential of {} MV since there is already {} processes of this potential running".format(potentia,max_nb_proc_per_pot))
      return False
    
    print("Starting potential of {} MV.".format(potentia))
    
    python_file = 'python_file_auto_' + str(a_random_num+iidx) + '.py'
  
    ## creating files
    batch_script_content = get_batch_script_content(python_file,buid_folder,nb_cpuu,potentia)
    wr = open(slurm_batch_file, 'w')  # opening and deleting
    wr.write(batch_script_content)
    wr.close()
    time.sleep(0.4)
    python_script_content = get_python_script_content(potentia,nb_shoo,nb_cpuu)
    wr = open(python_file, 'w')  # opening and deleting
    wr.write(python_script_content)
    wr.close()
    time.sleep(0.4)

    ## launching bash commands
    bashCommand = "sbatch " + slurm_batch_file
    process = sp.Popen(bashCommand.split(), stdout=sp.PIPE)
    output, error = process.communicate()

    print('output: {}'.format(output))
    print('Error: {}'.format(error))
    time.sleep(0.25)
    return True

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

if __name__ == "__main__":

  a_random_number = random.randint(1,100000000)

  minimum_record_per_potential = 500000

  process_ID_lst, potential_lst = make_running_id_and_potential_list() # get ID list and the corresponding potential list

  print("  ")
  print("process list ID:")
  print(process_ID_lst)
  print("corresponding potentials:")
  print(potential_lst)  
  print("  ")
  
  #potential_list =  [25,    -25,   50,   -50,  100, -100, 165, -165, 200, -200, 230, -230, 147.5, -147.5]
  #nb_shoot_list  =  [5000,  5000,  2500, 2500, 300,  300,  35,   35,  15,   15,  20,   20,  5000,   5000]
  #nb_cpu_list    =  [32,     32,   32,    32,   32,   32,   1,    1,   1,    1,   1,    1,     1,      1]
  
  # set of potentials and parameters for the runs on the cluster
  #potential_list =  [ -0.1,    50,   -50,      100,  -100,  150.,   -150,  200,   -200,  125,  -125,    250,   -250,  175.,  -175.,  225.,   -225.]
  #nb_shoot_list  =  [ 500,    200,   200,     5,     5,     2,      2,    2,      2,    2,     2,      1,      1,      1,     1,      1,       1.  ]
  #nb_cpu_list    =  [   32,    32,    32,      32,    32,    32,     32,   20,     20,   32,    32,     20,     20,     32,    32,     32,      32 ]

#-10.,10.,-30.,30.,-60.,60.,-80.,80.,-100.,100.,-150.,150.,-200.,200.,-300.,300.,-400.,400.
  #potential_list =  [   -3.,   3.,   -60.,   60.,   -100.,  100.,  -150.,  150.,  -200.,  200.,   -300.,  300.,  -250.,  250.,  325.,  -325., 360., -360.,   420.]
  #nb_shoot_list  =  [   2000,  2000, 2000,  2000,    2000,  2000,   2000,  2000,   2000,  2000,   2000,  2000,    2000,  2000,  2000,  2000 , 2000,  2000,   2000]
  #nb_cpu_list    =  [   32,    32,     32,    32,      32,    32,     32,    32,     32,    32,     32,    32,      32,  32  ,    32,    32 , 32,      32,   32]

  potential_list =  [   -3.,   3.,   -100.,  100.,  -150.,  150.,  -200.,  200.,  -300.,  300.,   -250.,  250.,  325.,  -325., 360., -360.,   420.]
  nb_shoot_list  =  [   2000,  2000,  2000,  2000,   2000,  2000,   2000,  2000,   2000,  2000,    2000,  2000,  2000,  2000 , 2000,  2000,   2000]
  nb_cpu_list    =  [   32,    32,      32,    32,     32,    32,     32,    32,     32,    32,      32,  32  ,    32,    32 , 32,      32,   32]


  #potential_list =  [   -3.,   3.,  -30.,  30.,  -60.,   60.,  -80.,   80.,  -100.,  100.,  -150.,  150.,  -200.,  200.,   -300.,  300.,  -250.,  250., 325., -325., 350., -350. ]
  #nb_shoot_list  =  [   3000,  3000,  3000, 3000,  3000,  3000,  3000,  3000,   3000,  3000,   3000,  3000,   3000, 3000,  3000,  3000,    3000,  3000, 3000, 3000,  3000,  3000 ]
  #nb_cpu_list    =  [   32,    32,    32,    32,   32,    32,    32,    32,     32,    32,     32,    32,     32,    32,     32,    32,      32,  32  ,   32,   32,    32,   32  ]
  
  statistics_dictionary = check_statistics_v2(potential_list)
  
  for key in statistics_dictionary:
      print("{} : {}".format(key,statistics_dictionary[key]))

  buid_folder = '/cluster/work/users/dsarria/G4_projects/BELOW_new_method/build/'

  slurm_batch_file = 'batch_file_auto.sh'

  for idx,potential in enumerate(potential_list):
    print("  ")
    key = str(int(float(potential)))

    if (statistics_dictionary[key]<minimum_record_per_potential):
        print("Stats of {} for {} MV is less than the requirement of {}.".format(statistics_dictionary[key],potential,minimum_record_per_potential))
        done_or_not = launch_job(potential, nb_shoot_list[idx], nb_cpu_list[idx], a_random_number, idx, process_ID_lst, potential_lst)
    else:
        print("Stats of {} for {} MV is OK.".format(statistics_dictionary[key],potential))
    
  print("  ")




