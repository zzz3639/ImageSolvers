delete ./MexFile/*
mex Ccode/EMFixSmoothMex.c -outdir ./MexFile
mex Ccode/EMSparseSmoothMex.c -outdir ./MexFile
mex Ccode/MergeMex.c -outdir ./MexFile