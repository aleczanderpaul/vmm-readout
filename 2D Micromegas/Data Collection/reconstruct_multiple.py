#!/usr/bin/python

import os
import subprocess
import re
import sys
import glob

pcapngFolder = "/home/r2d2/vmm-sdat/data-taking/VMMDataTestpcapngFilesMay7"
pcapngFiles = glob.glob(os.path.join(pcapngFolder, "*.pcapng"))

for filePath in pcapngFiles:
    try:
	    args = ['/home/r2d2/vmm-sdat/build/convertFile',
	    '-f', f'{filePath}',
	    '-geo', '/home/r2d2/vmm-sdat/reconstruction/small_detector_updated.json',
	    '-bc', '40', '-tac', '60',
	    '-th', '0',
	    '-cs', '2', '-ccs', '4', '-mst', '1',
	    '-dt','200','-spc','1500','-dp', '200',
	    '-crl','0',
	    '-cru','1000',
	    '-save', '[[1],[1],[1]]',
	    '-info', 'test',
	    '-df', 'SRS', '-cal','/home/r2d2/vmmsc/calibs/20250505/9mVfC_200ns_60ns_negative_321DAC/vmm_calibration_time_FEC6_VMM0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_123448.json']
	    subprocess.call(args)   

    except OSError:
	    pass
