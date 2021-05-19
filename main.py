
import splas
import numpy as np
import splasIO as iojson
import sys
import os
from sys import argv
from subprocess import Popen
from mpi4py import MPI
import time

################################################################################################################################################################################################################

def ensemble_ns():

    file1 = "./input/ns0.01.json"
    file2 = "./input/ns0.1.json"
    file3 = "./input/ns0.5.json"
    filenew1 = "./input/ns2.json"
    filenew2 = "./input/ns3.5.json"
    filenew3 = "./input/ns7.json"
    file4 = "./input/ns5.json"
    filedef = './input/input0.json'
    file5 = "./input/ns10.json"
    file6 = "./input/ns20.json"
    ns = [file1,file2,file3,filenew1,filenew2,filenew3,file4,filedef,file5,file6]
    return ns

################################################################################################################################################################################################################

def ensemble_xi():

    file1 = "./input/zxi005.json"
    file2 = "./input/zxi015.json"
    filedef = "./input/input0.json"
    file3 = "./input/zxi01.json"
    file4 = "./input/zxi025.json"
    xi = [file1,file2,filedef,file3,file4]
    return xi

################################################################################################################################################################################################################

def ensemble_dx():

    file1 = "./input/dx13.json"
    file2 = "./input/dx015.json"
    file3 = "./input/dx016.json"
    file4 = "./input/dx0177.json"
    filedef = './input/input0.json'
    file5 = "./input/dx0354.json"
    file6 = "./input/dx05.json"
    file7 = "./input/dx1.json"
    file8 = "./input/dx02.json"
    dx = [file1,file2,file3,file4,filedef,file5,file6,file7,file8]
    return dx

################################################################################################################################################################################################################
################################################################################################################################################################################################################

def ensemble_ni():

    file1 = "./inputIC/n1.json"
    file2 = "./inputIC/n2.json"
    file3 = "./inputIC/n3.json"
    file4 = "./inputIC/n4.json"
    file5 = "./inputIC/n5.json"
    file6 = "./inputIC/n6.json"
    filedef = "./inputIC/input0.json"
    ni = [file1,file2,file3,file4,file5,file6,filedef]
    return ni

################################################################################################################################################################################################################

def ensemble_pi():

    file1 = "./inputIC/p1.json"
    file2 = "./inputIC/p2.json"
    file3 = "./inputIC/p3.json"
    file4 = "./inputIC/p4.json"
    file5 = "./inputIC/p5.json"
    file6 = "./inputIC/p6.json"
    filedef = "./inputIC/input0.json"
    pi = [file1,file2,file3,file4,file5,file6,filedef]
    return pi

################################################################################################################################################################################################################

def ensemble_zi():

    file1 = "./inputIC/z1.json"
    file2 = "./inputIC/z2.json"
    file3 = "./inputIC/z3.json"
    file4 = "./inputIC/z4.json"
    file5 = "./inputIC/z5.json"
    file6 = "./inputIC/z6.json"
    filedef = "./inputIC/input0.json"
    zi = [file1,file2,file3,file4,file5,file6,filedef]
    return zi

#############################################################################################################################################################################################
#############################################################################################################################################################################################

def ensemble_qi():

    file1 = "./sink0.json"
    file2 = "./sink1.json"
    file3 = "./sink2.json"
    file4 = "./sink3.json"
    file5 = "./sink4.json"
    filedef = "./default.json"
    file6 = "./sink5.json"
    qi = [file1,file2,file3,file4,file5,filedef,file6]
    return qi

#############################################################################################################################################################################################
#############################################################################################################################################################################################

def ensemble_pmort():

    file1 = "./input/pmort004.json"
    file2 = "./input/pmort005.json"
    filedef = "./input/input0.json"
    file3 = "./input/pmort02.json"
    file4 = "./input/pmort04.json"
    file5 = "./input/pmort06.json"
    file6 = "./input/pmort075.json"
    file7 = "./input/psink0.json"
    file8 = "./input/psink0.001.json"
    file9 = "./input/psink0.01.json"
    pmort = [file1,file2,filedef,file3,file4,file5,file6,file7,file8,file9]
    return pmort

#############################################################################################################################################################################################
#############################################################################################################################################################################################

def ensemble_zpar_test():

    filedef = "./input_test_1.5/input0.json"
    file1 = "./input_test_1.5/eps06.json"
    file2 = "./input_test_1.5/feg06.json"
    file7 = "./input_test_1.5/kz1.json"
    file8 = "./input_test_1.5/kz2.json"
    file9 = "./input_test_1.5/kz4.json"
    file10 = "./input_test_1.5/kz5.json"
    otest = [filedef,file1,file2,file7,file8,file9,file10]
    return otest

#############################################################################################################################################################################################
#############################################################################################################################################################################################

def ensemble_double():

    #file1 ="./sensitivity_json/sensitivity_diff1.json"
    #file2 ="./sensitivity_json/sensitivity_diff2.json"
    #file3 ="./sensitivity_json/sensitivity_wp1.json"
    #file4 ="./sensitivity_json/sensitivity_wp2.json"
    #double = [file1,file2,file3,file4]
    file1 = "./input/default4.json"
    file2 = "./input/gumble4.json"
    simu = splas.splas(file1)
    simu.get_input()
    simu.write()
    simu.sim()
    simu.store()
    print("fine sim1")
    simu2 = splas.splas(file2)
    simu2.get_input()
    simu2.write()
    simu2.sim()
    simu2.store()
    print("fine sim2")


#############################################################################################################################################################################################
#############################################################################################################################################################################################

def one_sim():

    filexx = "./input_ps1/input6.json"
    #filexx = "./input_closed_world_tv/input1.json"
    #filexx = "./input_closed_world_tv/input2.json"
    #filexx = "./input_closed_world_tv/input3.json"
    #filexx = "./input_closed_world_tv/input4.json"
    #filexx = "./input_closed_world_tv/input5.json"
    #filexx = "./input_closed_world_tv/input6.json"
    print("Using file {}".format(filexx))
    simu = splas.splas(filexx)
    simu.get_input()
    simu.write()
    simu.sim()
    simu.store()
    #print('phyto 1'+ str(simu.phyto[1,1,1]))
    #print('phyto 2'+str(simu.phyto[1,2,1]))
    #print('det 1'+str(simu.detritus[1,1]))
    #print('N_1'+str(simu.Nitrogen[1,20]))
    #print('det_100'+str(simu.detritus[100,20]))
    #simu.get_simulation_output()

################################################################################################################################################################################################################

def master(comm,worklist):
    print("START MASTER ",comm.Get_rank())
    st = MPI.Status()
    slave_num = 0
    for simfile in worklist:
        slave_num = comm.recv(source=MPI.ANY_SOURCE,tag=0, status=st)
        comm.send(simfile, dest=slave_num, tag=1)
        print("SEND TO:",slave_num," FILE = ",simfile)

    total_num_slaves = comm.Get_size()-1
    while(total_num_slaves>0):
        slave_num = comm.recv(source=MPI.ANY_SOURCE,tag=0, status=st)
        comm.send("_STOP_", dest=slave_num, tag=1)
        total_num_slaves = total_num_slaves - 1

    print(comm_rank," MASTER END COMPUTATION ")


def slave(comm):
    print("START SLAVE ",comm.Get_rank())
    active_slave = True
    comm_rank = comm.Get_rank()
    st = MPI.Status()

    while(active_slave):
        dimension = 0
        comm.send(comm_rank, dest=0, tag=0)
        filedef = comm.recv(source=0,tag=1, status=st)

        if(filedef =="_STOP_"):
            print(comm_rank," KILLED ")
            active_slave = False
            continue

        #recv
        t1=time.time()
        default = splas.splas(filedef)
        default.get_input()
        default.write()
        default.sim()
        default.store()

        #write files for results
        default.get_simulation_output()

        t2=time.time()
        print(comm_rank,filedef," time = ",t2-t1)
        del(default)


#MAIN
if __name__ == '__main__':

    worklist = []
    comm = MPI.COMM_WORLD
    comm_rank = comm.Get_rank()

    if(comm_rank == 0):

        try:
            script, simflag = argv
            if simflag == '1' or simflag == 'One' or simflag == 'one':
                print("One single simulation")
                one_sim()
            elif simflag == 'ns':
                print("Ensemble simulation varying Nutrient supply")
                worklist = ensemble_ns()
            elif simflag == 'xi':
                print("Ensemble simulation varying Zooplankton mortality")
                worklist = ensemble_xi()
            elif simflag == 'dx':
                print("Ensemble simulation varying Zooplankton prey selectivity")
                worklist = ensemble_dx()
            elif simflag == 'ni':
                print("Ensemble simulation varying initial Nitrogen concentration")
                worklist = ensemble_ni()
            elif simflag == 'pi':
                print("Ensemble simulation varying initial Phytoplankton abundance")
                worklist = ensemble_pi()
            elif simflag == 'zi':
                print("Ensemble simulation varying initial Zooplankton abundance")
                worklist = ensemble_zi()
            elif simflag == 'qi':
                print("Ensemble simulation varying initial Phytoplankton Cell Quota (#individuals)")
                worklist = ensemble_qi()
            elif simflag == 'pmort':
                print("Ensemble simulation varying Phytoplankton mortality")
                worklist = ensemble_pmort()
            elif simflag == 'ztest':
                print("Ensemble simulation varying some zooplankton parameters (k,eps,feg)")
                worklist = ensemble_zpar_test()
            elif simflag == 'double':
                print("Double simulation")
                worklist = ensemble_double()
            else:
                print("Simulation not specified correctly")
        except ValueError:
            print("""Please specify which simulation you'd like to perform.
                     You must run 'pythonV main.py flag'
                     where V is python version used and 'flag' is:

                                    - '1' for one single simulation
                                    - 'ns' for an ensemble simulation varying nutrient load
                                    - 'xi' for an ensemble simulation varying all zooplankton mortality
                                    - 'dx' for an ensemble simulation varying zooplankton prey selectivity
                                    - 'ni' for an ensemble simulation varying initial N concentration
                                    - 'pi' for an ensemble simulation varying initial P abundance
                                    - 'zi' for an ensemble simulation varying initial Z abundance
                                    - 'qi' for an ensemble simulation varying initial P cell quota (#ind)
                                    - 'pmort' for ensemble simulation varying phytoplankton mortality
                                    - 'ztest' for ensemble simulation varying some zooplankton parameters
                                    - 'double' for a double simulation in parallel

        ---------------------------------------------------------------------------------------------
        ---------------------------------------------------------------------------------------------
            IF YOU GET THIS ERROR IN THE MIDDLE OF THE EXECUTION IT MEANS THAT THERE IS A ValueError,
            SO YOU SHOULD INVESTIGATE FURTHER.
            IF NOT, PLEASE SIMPLY INSERT THE DESIRED SIMULATION.

            """)

        master(comm,worklist)

    else:

        slave(comm)

    #comm.Barrier()
