%
% casenum = 1 call: yfun.m    [obtain realizations of the uncertain subsystem, y(a,e,t)]
% casenum = 2 call: zfun.m    [obtain realizations of the integrated system, z1(a,e,theta,t), z2(a,e,theta,t)]
% casenum = 3 call: gfun.m    [compute the requirements vector g(a,e,theta)]
% casenum = 4 call: yfun.m + zfun.m  [obtain realizations uncertain subsystem and integrated system using a single system command]
% casenum = 5 call: yfun.m + gfun.m  [obtain realizations uncertain subsystem and requirements vector using a single system command]
% casenum = 6 call: zfun.m + gfun.m  [obtain realizations integrated system and requirements vector using a single system command]
% casenum = 7 call: yfun + zfun.m + gfun.m  [obtain realizations of y,z,g using a single system command]
%


afile='aleatory.dat';
efile='epistemic.dat';
dfile='design.dat';

na=5;       % Number of aleatory variables  (cannot be changed)
ne=4;       % Number of epistemic variables (cannot be changed)
nd=9;       % Number of design variables    (cannot be changed)
%
% Nr was created for use with the version based upon the Matlab Runtime Library to decrease the impact of the
% lenghty start-up delay. The use of "nr" for the pure C standalone version not needed, but the functionality has 
% been retained. We suggest using nr=1;
%
nr=20;      % Number of random realizations of aleatory, epistemic, and design parameters (user defined)
nr=1;       % Number of random realizations of aleatory, epistemic, and design parameters (user defined)

rng default                 % Set default random number generator

%
% As an example only -  just assign random numbers.....
%
a=rand(nr*na,1);
e=rand(nr*ne,1);
theta=rand(1,nr*nd);


fileAID = fopen(afile,'w');
fileEID = fopen(efile,'w');
fileDID = fopen(dfile,'w');

fprintf(fileAID,'%.20f\n',a);
fprintf(fileEID,'%.20f\n',e);
fprintf(fileDID,'%.20f\n',theta);

fclose(fileAID);
fclose(fileEID);
fclose(fileDID);

