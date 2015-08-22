function [ pic ] = PostRun( pic0, ProcessPara )
%POSTRUN Summary of this function goes here
% Usage: [ pic ] = PostRun( pic0, ProcessPara )
%   pic0: Molecule table to be processed
%   ProcessPara: [MergeDist, MolZero, FoldRemain]
%      MergeDist: Molecules within MergeDist are treated as one
%      MolZero: Molecules below this intensity are ignored.
%      FoldRemain: [0-1], Molecules are solted by intensity, only molecules
%                  contribute to top FoldRemain are preserved.

MergeDist=ProcessPara(1);
MolZero=ProcessPara(2);
FoldRemain=ProcessPara(3);
pic=Merge(pic0, MergeDist, MolZero);
pic=removezeros(pic,FoldRemain);

end
