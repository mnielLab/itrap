#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
import shutil
import subprocess
import tempfile
import threading
from contextlib import contextmanager

seq2score="/home/tuba/herpov/tcr-pmhc-sc-project/tools/seq2score_db_kernel"
BLF="/home/tuba/herpov/tcr-pmhc-sc-project/tools/BLOSUM50"
QIJ="/home/tuba/herpov/tcr-pmhc-sc-project/tools/blosum62.qij"

F1 = "/home/tuba/herpov/tcr-pmhc-sc-project/notebooks/lstS.txt"
F2 = "/home/tuba/herpov/tcr-pmhc-sc-project/notebooks/lstM.txt"

lstS = pd.DataFrame(['CAVRSAYSGAG'])
lstM = pd.DataFrame(['CAARLIQGAQKLVF', 'CAGPSYNTDKLIF', 'CAMPNSGGYQKVTF', 'CAMNRDDKIIF', 'CAVRSAYSGAGSYQLTF'])


cmd = [seq2score, '-blf', BLF, '-blqij', QIJ, '-pa', F1, F2] #

def run_cmd(cmd, input_string=''):
    """Run the cmd with input_string as stdin and return output."""
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                                             stderr=subprocess.PIPE, universal_newlines=True, close_fds=True)
    out, stderr = p.communicate(input=input_string)
    if p.returncode:
            raise Exception('Cmd {} failed: {}'.format(cmd[0], stderr))
    return out

#print(run_cmd(cmd))


@contextmanager
def named_pipes(count):
    dirname = tempfile.mkdtemp()
    try:
        paths = []
        for i in range(count):
            paths.append(os.path.join(dirname, 'named_pipe' + str(i)))
            os.mkfifo(paths[-1])
        yield paths
    finally:
        shutil.rmtree(dirname)

def write_command_input(lst, path):
    lst.to_csv(path, header=False,index=False)
    #np.savetxt(path, lst, fmt='%s')

def run_seq2score():
	with named_pipes(2) as paths:
	    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
	                                         stderr=subprocess.PIPE, universal_newlines=True, close_fds=True) #
	    with p.stdout:
	        for lst, path in zip([lstS, lstM], paths):
	            t = threading.Thread(target=write_command_input, args=[lst, path]) 
	            t.daemon = True
	            t.start()
	        result = pd.read_csv(p.stdout, sep=" ", names=['seq1', 'seq2', 'similarity'], usecols=[1,2,3], comment='#')   
	p.wait()
	return result



#print(run_seq2score())

def run_cmd(cmd, input_string=''):
    """Run the cmd with input_string as stdin and return output."""
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                                             stderr=subprocess.PIPE, universal_newlines=True, close_fds=True)
    result = pd.read_csv(p.stdout, sep=" ", names=['seq1', 'seq2', 'similarity'], usecols=[1,2,3], comment='#')
    p.wait()
    return result

print(run_cmd(cmd))

