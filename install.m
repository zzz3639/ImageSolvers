delete ./MexFile/*
mex Ccode/EMFixSmoothMex.c -outdir ./MexFile
mex Ccode/EMSparseSmoothMex.c -outdir ./MexFile
mex Ccode/MergeMex.c -outdir ./MexFile

mex Ccode/EMMain.c Ccode/EMfunctions.c -DPosFix=false -outdir ./MexFile -output EMSparseMex
mex Ccode/EMMain.c Ccode/EMfunctions.c -DPosFix=true -outdir ./MexFile -output EMFixMex
mex Ccode/EMRobustMain.c Ccode/EMfunctions.c -outdir ./MexFile -output EMRobustMex
%mex Ccode/try.c Ccode/EMfunctions.c -outdir ./MexFile -output tryMex