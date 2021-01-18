# run using: python test.py
#
#  This script presents an example of calling the standalone C version
#  of the 2019 NASA UQ Challenge Problem from within python.
#

import numpy as np
import matplotlib.pyplot as plt
import os
import csv


#
# ---------------------------------------------------------------------------------------------------
#  Run UQ Challenge Probelm
#
#  The variable casenum is defined in datafile: "casefile.dat"
#
# casenum = 1 call: yfun.m    [obtain realizations of the uncertain subsystem, y(a,e,t)]
# casenum = 2 call: zfun.m    [obtain realizations of the integrated system, z1(a,e,theta,t), z2(a,e,theta,t)]
# casenum = 3 call: gfun.m    [compute the requirements vector g(a,e,theta)]
# casenum = 4 call: yfun.m + zfun.m  [obtain realizations uncertain subsystem and integrated system using a single system command]
# casenum = 5 call: yfun.m + gfun.m  [obtain realizations uncertain subsystem and requirements vector using a single system command]
# casenum = 6 call: zfun.m + gfun.m  [obtain realizations integrated system and requirements vector using a single system command]
# casenum = 7 call: yfun + zfun.m + gfun.m  [obtain realizations of y,z,g using a single system command]
#
# The following command runs the standalone Matlab executable for casenum = 7 
#
cmd = ' uqsim 7  "aleatory.dat" "epistemic.dat" "design.dat"'
os.system(cmd)

#
# Read response data from datafiles
#

t = np.genfromtxt('tout')
yout = np.genfromtxt('yout')
yout = yout.transpose()
plt.figure()
plt.plot(t,yout)
plt.ylabel('Predicted y(t)')
plt.xlabel('Time, [s]')
plt.show(block = False)

z1out = np.genfromtxt('z1out')
z2out = np.genfromtxt('z2out')
z1out = z1out.transpose()
z2out = z2out.transpose()

plt.figure()
plt.plot(t,z1out)
plt.ylabel('Predicted z_1(t)')
plt.xlabel('Time, [s]')
plt.show(block = False)

plt.figure()
plt.plot(t,z2out)
plt.ylabel('Predicted z_2(t)')
plt.xlabel('Time, [s]')
plt.show(block = False)

plt.show()


