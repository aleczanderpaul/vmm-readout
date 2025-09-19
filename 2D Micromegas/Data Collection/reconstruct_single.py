#!/usr/bin/python

import os
import subprocess
import re
import sys

try:
	args = ['/home/tpcmaster/vmm-sdat/build/convertFile',
	'-f', '/home/tpcmaster/runs/AUG23/DLC/Fe55/700Vmesh.pcapng',
	'-geo', '/home/tpcmaster/runs/AUG23/small_detector_geo.json',
	'-bc', '40', '-tac', '60',
	'-th', '0',
	'-cs', '2', '-ccs', '4', '-mst', '1',
	'-dt','200','-spc','1500','-dp', '200',
	'-crl','0',
	'-cru','1000',
	'-save', '[[1],[1],[1]]',
	'-info', 'test',
	'-df', 'SRS', '-cal','/home/tpcmaster/runs/AUG23/calibs/ADC_Time_DLC/AUG23_DLC_calib.json']
	subprocess.call(args)

except OSError:
	pass
