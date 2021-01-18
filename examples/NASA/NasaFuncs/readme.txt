Usage:

uqsim casenum  aleatory.dat epistemic.dat design.dat   (casenum is an integer, and *.dat are floating point datafiles)

#
# casenum = 1 call: yfun.m    [obtain realizations of the uncertain subsystem, y(a,e,t)]
# casenum = 2 call: zfun.m    [obtain realizations of the integrated system, z1(a,e,theta,t), z2(a,e,theta,t)]
# casenum = 3 call: gfun.m    [compute the requirements vector g(a,e,theta)]
# casenum = 4 call: yfun.m + zfun.m  [obtain realizations uncertain subsystem and integrated system using a single system command]
# casenum = 5 call: yfun.m + gfun.m  [obtain realizations uncertain subsystem and requirements vector using a single system command]
# casenum = 6 call: zfun.m + gfun.m  [obtain realizations integrated system and requirements vector using a single system command]
# casenum = 7 call: yfun + zfun.m + gfun.m  [obtain realizations of y,z,g using a single system command]
#

Example: compute yfun, zfun and gfun in a single function call

uqsim 7 aleatory.dat epistemic.dat design.dat

Additional Resources:

See test.py for an example of using python with uqsim

Also see: make_data.m for an example of creating the required input files (aleatory.dat, epistemic.dat, and design.dat). 
The file make_data.m is a Matlab script, so users will need to develop an equivalent for their chosen coding environment.

