function [FrameInd CAP_Ind]=Comp_FrameInd(Params,TH,Cap_par)
 
% This function extracts Time Frames and CAPs indices from your params. 
% Needed for the stats (Amico et al., PloS one 2014).
%  
% usage:
% [FrameInd CAP_Ind] = Comp_FrameInd(Params,TH,Cap_par)   
%
% Input:
% Params - Struct array of parameters computed using the Comp_Params
% function
% TH - the threshold for the time frames collection (e.g. in the paper, 1SD
% is ~= 15%, so in that case TH=15). It can be a single number or an array. 
% Cap_par - Struct array containing information about each CAP computed
% using the Comp_CAP function.
%
% Output: 
% FrameInd - Struct array of the temporal indices for the time Frames extracted
% CAP_Ind - Struct array of CAP labels for every time Frame collected.
% Ref: Amico, Enrico, et al. "Posterior Cingulate Cortex-Related 
% Co-Activation Patterns: A Resting State fMRI Study in Propofol-Induced 
% Loss of Consciousness." PloS one 9.6 (2014): e100012.
% 
% Written by EA Sep 2014.
%__________________________________________________________________________
% License
%
% This file is part of CAPsTOOLBOX.  It is Copyright (C) Enrico Amico, 
% 
% CAPsTOOLBOX is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% CAPsTOOLBOX is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% For a copy of the GNU General Public License, see <http://www.gnu.org/licenses/>.
%__________________________________________________________________________

Perc = Params.Rate1;
Index = Params.Ind1;
CAP_Ind = cell(size(Index,1),length(TH));
FrameInd = cell(size(Index,1),length(TH));
for pp = 1:length(TH) 
    thr = TH(pp); %%% percentage of data to consider
    pick=1;
    CInd = Cap_par.ClusterInd;
    for i=1:size(Index,1)         
        Ix = find(Perc(i,:)>=(100-thr)); 
        IndFrames = Index{i,Ix(end)}; %% time points of the frames 
        CAP_Ind{i,pp} = CInd(pick:pick+length(IndFrames)-1); %% CAP index
        FrameInd{i,pp} = IndFrames; 
        pick = pick+length(IndFrames); 
    end
    clear Ix;
end

return;
 
