function [Pic, No, Mv] = RunSolverTunning(Img, NumMol, Lambda, OptPara, Tolerance1, Tolerance2, ProcessPara1, ProcessPara2)
%Modified by ZHANG Haowen in 2015.08.29
%
% Usage: [Pic,No,Mv]=RunSolverTunning(Img,NumMol,Lambda,OptPara,Tolerance1,Tolerance2,ProcessPara1,ProcessPara2)
% Input variables:
%   Img: m by n matrix, raw image
%   NumMol: number of initial molecules, before merge. Integer
%   Lambda: Lambda for the optimization.
%   OptPara: [Sigma, BoundarySize, PSFdecay];
%   Tolerance1: [MaxIte0, MaxIteFast, Izero1, Pzero1] for FastSolver
%      MaxIte0: Maximal iteration number of the first SparseGMM run of FastSolver
%      MaxIteFast: Maximal iteration number of the second SparseGMM run of FastSolver
%      Izero1: minimal Intensity change of FastSolver
%      Pzero1: minimal Position change of FastSolver
%   Tolerance2: [MaxIteTunning, Izero2, Pzero2] for RunSolverFix
%      MaxIteTunning: Maximal iteration number of RunSolverFix
%      Izero2: minimal Intensity change of RunSolverFix
%      Pzero2: minimal Position change of RunSolverFix
%   ProcessPara1: [MergeDist1, MolZero1] for FastSolver
%      MergeDist1: Molecules within MergeDist are treated as one
%      MolZero1: Molecules below this intensity are ignored.
%   ProcessPara2: [MergeDist2, MolZero2] for PostRun
%      MergeDist2: Molecules within MergeDist are treated as one
%      MolZero2: Molecules below this intensity are ignored.
%      FoldRemain: [0-1], Molecules are solted by intensity, only molecules
%                  contribute to top FoldRemain are preserved.


% Run FastSolver
[pic, no, Mv] = FastSolver(Img, NumMol, Lambda, OptPara, Tolerance1, ProcessPara1);
pic=PostRun(pic,ProcessPara2);

% Run Tunning
mv0.pic=pic;
mv0.no=no;
[Pic,No]=RunSolverFix( Img, OptPara, Tolerance2, mv0 );

end

