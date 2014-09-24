function data = Comp_LiuMask(data,dim) 

% This function improves the signal to noise ratio over the time frames
% collected.
%  
% usage:
% data = Comp_LiuMask(data,dim);
%
% Input:
% data - NxM matrix for N time points (e.g. functional volumes) in M voxels;
% dim - the x, y and z dimensions of the volume (e.g. dim=[61 71 63])
%
% Output: 
% data - the data matrix after masking
%
% Ref: Amico, Enrico, et al. "Posterior Cingulate Cortex-Related 
% Co-Activation Patterns: A Resting State fMRI Study in Propofol-Induced 
% Loss of Consciousness." PloS one 9.6 (2014): e100012.
% 
% Liu, Xiao, and Jeff H. Duyn. "Time-varying functional network information 
% extracted from brief instances of spontaneous brain activity." Proceedings 
% of the National Academy of Sciences 110.11 (2013): 4392-4397.
% 
% Written by EA Sep 2014
%__________________________________________________________________________
%License
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


Nvol = size(data,1);
p1 = 0.05; %% threshold
Nvox = size(data,2);
for i = 1:Nvol
    datavec = data(i,:);
    t1 = p1 * max(datavec); %% 10% high value
    t2 = p1 * min(datavec); %% 10% low value
    datavec(datavec<t1 & datavec>0) = 0;
    datavec(datavec>t2 & datavec<0) = 0;
    tmp = reshape(datavec,dim); 
    L = bwlabeln(tmp,6);
    tmp(L==0) = 0;
    datavec = reshape(tmp,[1 Nvox]);
    data(i,:) = datavec;
end

return;
