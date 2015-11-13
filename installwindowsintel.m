delete ./MexFile/*
mex Ccode/EMFixSmoothMex.c -outdir ./MexFile COMPFLAGS="$COMPFLAGS /Qstd=c99"
mex Ccode/EMSparseSmoothMex.c -outdir ./MexFile COMPFLAGS="$COMPFLAGS /Qstd=c99"
mex Ccode/MergeMex.c -outdir ./MexFile COMPFLAGS="$COMPFLAGS /Qstd=c99"

mex Ccode/EMMain.c Ccode/EMfunctions.c -DPosFix=false -outdir ./MexFile -output EMSparseMex COMPFLAGS="$COMPFLAGS /Qstd=c99"
mex Ccode/EMMain.c Ccode/EMfunctions.c -DPosFix=true -outdir ./MexFile -output EMFixMex COMPFLAGS="$COMPFLAGS /Qstd=c99"
mex Ccode/EMRobustMain.c Ccode/EMfunctions.c -outdir ./MexFile -output EMRobustMex COMPFLAGS="$COMPFLAGS /Qstd=c99"
mex Ccode/EMMainUneven.c Ccode/EMfunctions.c -DPosFix=false -outdir ./MexFile -output EMSparseUnevenMex COMPFLAGS="$COMPFLAGS /Qstd=c99"
mex Ccode/EMMainUneven.c Ccode/EMfunctions.c -DPosFix=true -outdir ./MexFile -output EMSparseUnevenFixMex COMPFLAGS="$COMPFLAGS /Qstd=c99"