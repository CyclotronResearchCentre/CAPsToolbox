function [Params] = Comp_Params(TS,Data,CMap,brind,iseed) 

% This function computes the params needed for the CAP computation.
%  
% usage:
% [Params] = Comp_Params(TS,Data,CMap,brind,iseed)
%
% Input:
% TS - cell array (number of subjects,number of seeds). Each cell contain the seed Time 
% course Mx1 array with M=number of time points (e.g fMRI volumes).
% Data - cell array (1,number of subjects). Each cell should contain an 
% NxM matrix for N time points (e.g. functional volumes) in M voxels;
% CMap - cell array (number of subjects, number of seeds). Each cell 
% contain the seed-voxel correlation map as a 1xM array with M=number of 
% voxels.It can be saved into an .img/.nii file using the spm_write 
% function (see spm_write help for info). 
% brind - brain mask, organized as 1xM, with M=number of brain voxels.
% iseed - seed index: it specifies which seed to consider (if there is just
% one seed, then iseed=1).
%
% Output: 
% Params= Struct array of parameters needed for CAPs computation.
%
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

Params ={};

for i=1:length(Data)
    fprintf('\n Processing subject %d \n',i);  
    CorrMap = CMap{i,iseed}; 
    New_data = Data{i};
    New_data(isnan(New_data))=0;
    Tseries = TS{i,iseed};
    TotN = max(size(Tseries));
    [M indM] = max(Tseries);
    [m indm] = min(Tseries);
    Start = M;
    Step = (M-m)/50; 
    count = 1;
    fprintf('\n Threshold up %d \n',i);  
    for thr = Start:-Step:m
        fprintf('.');
        index1 = find(Tseries >= thr);
        IndStruct{i,count} = index1; 
        ActiMap1 =zeros(1,size(New_data,2));
        tmp = mean(New_data(index1,:),1); %% something like 1 x nvoxels
        ActiMap1(brind) = tmp(brind);
        rate1(i,count) = (1-(length(index1)/TotN))*100;
        SpatCorr1(i,count) = corr(CorrMap',ActiMap1'); %%% needs transpose for this kind of arrays
        count = count+1;
    end
    
     
      

end

Params.Ind1 = IndStruct;
Params.Rate1 = rate1;
Params.SpatCorr1 = SpatCorr1;
SpatCorr1M = mean(SpatCorr1,1);
rate1M = mean(rate1,1);
figure;
plot(rate1M,SpatCorr1M);
set(gca,'xdir','reverse');

return;
