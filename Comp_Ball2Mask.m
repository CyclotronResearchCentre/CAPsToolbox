function [cor_s] = Comp_Ball2Mask(ABrainSize, AVoxelSize, AROICenter, AROIRadius, Header)

% This function creates spherical seeds from your image
%  
% usage:
% CorrMap = Comp_CMap(data,V,brind,seed_mni,seed_name,seed_radius);
%
% Input:
% ABrainSize - the x, y and z dimensions of the volume (e.g. dim=[61 71 63])
% AVoxelSize - the x, y and z voxel size
% AROICenter - 3D coordinates (in mm) of your seed centroid (e.g seed_mni= [0 53 26])
% AROIRadius = your seed radius (in mm, e.g. seed_radius=6). 
% Header - spm vector of structures containing image volume information. Needed 
% as reference for the voxel to world transformation (see spm_vol for info)
% 
% Output: 
% cor_s = Matrix which contains the 3D coordinates of the voxels included
% in the spherical seed.
%
% Ref: Amico, Enrico, et al. "Posterior Cingulate Cortex-Related 
% Co-Activation Patterns: A Resting State fMRI Study in Propofol-Induced 
% Loss of Consciousness." PloS one 9.6 (2014): e100012.
% 
% Written by Guarong Wu 2010, updated by EA Sep 2014
%__________________________________________________________________________
%License
%
% This file is part of CAPsTOOLBOX. It is Copyright (C) Enrico Amico. 
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



AROICenter=reshape(AROICenter, 1,length(AROICenter));

if isfield(Header,'mat')
    AROICenter=round(inv(Header.mat)*[AROICenter,1]'); % From  world coordinates to voxel space Vspace = A^-1 Realspace
    AROICenter=AROICenter(1:3);
    AROICenter=reshape(AROICenter, 1,length(AROICenter));
else
    AROICenter =round(-1*AROICenter./AVoxelSize) +AOrigin;
end

radiusX =round(AROIRadius /AVoxelSize(1));% From mm to voxel size
if (AROICenter(1)-radiusX)>=1 && (AROICenter(1)+radiusX)<=ABrainSize(1) 
    rangeX	=(AROICenter(1)-radiusX):(AROICenter(1)+radiusX);
elseif (AROICenter(1)-radiusX)<1 && (AROICenter(1)+radiusX)<=ABrainSize(1)
    rangeX	=1:(AROICenter(1)+radiusX);
elseif (AROICenter(1)-radiusX)>=1 && (AROICenter(1)+radiusX)>ABrainSize(1)
    rangeX	=(AROICenter(1)-radiusX):ABrainSize(1);
else
    rangeX =1:ABrainSize(1);
end

radiusY =round(AROIRadius /AVoxelSize(2));
if (AROICenter(2)-radiusY)>=1 && (AROICenter(2)+radiusY)<=ABrainSize(2)
    rangeY	=(AROICenter(2)-radiusY):(AROICenter(2)+radiusY);
elseif (AROICenter(2)-radiusY)<1 && (AROICenter(2)+radiusY)<=ABrainSize(2)
    rangeY	=1:(AROICenter(2)+radiusY);
elseif (AROICenter(2)-radiusY)>=1 && (AROICenter(2)+radiusY)>ABrainSize(2)
    rangeY	=(AROICenter(2)-radiusY):ABrainSize(2);
else
    rangeY =1:ABrainSize(2);
end

radiusZ =round(AROIRadius /AVoxelSize(3));
if (AROICenter(3)-radiusZ)>=1 && (AROICenter(3)+radiusZ)<=ABrainSize(3)
    rangeZ	=(AROICenter(3)-radiusZ):(AROICenter(3)+radiusZ);
elseif (AROICenter(3)-radiusZ)<1 && (AROICenter(3)+radiusZ)<=ABrainSize(3)
    rangeZ	=1:(AROICenter(3)+radiusZ);
elseif (AROICenter(3)-radiusZ)>=1 && (AROICenter(3)+radiusZ)>ABrainSize(3)
    rangeZ	=(AROICenter(3)-radiusZ):ABrainSize(3);
else
    rangeZ =1:ABrainSize(3);
end
cor_s = [];
for x=rangeX, for y=rangeY, for z=rangeZ,
    %Ball Definition, Computing within a cubic to minimize the time to be consumed
    
    if norm(([x, y, z] -AROICenter).*AVoxelSize)<=AROIRadius, % back in mm world space for the mni coordinates

        cor_s = [cor_s; [x, y, z]];
    end
end; end; end;
