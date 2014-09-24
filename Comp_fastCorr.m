function R = Comp_fastCorr(X, Y)

% This function computes Pearson Correlation between two vectors
%  
% usage:
% R = Comp_fastCorr(X, Y,flag);
%
% Input:
% X, Y - the input vectors
% 
% Output: 
% R = The Correlation matrix.
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

[numSamp1 ~] = size(X); % n*p1,n is the number of samples and p1 is the number
X = (X - repmat(mean(X), numSamp1, 1))./repmat(std(X, 0, 1), numSamp1, 1);
Y = (Y - repmat(mean(Y), numSamp1, 1))./repmat(std(Y, 0, 1), numSamp1, 1);
R = X' * Y / (numSamp1 - 1);
R(isnan(R))=0; %%% when std 0 can happen to have Nan, for convention put it to zero.
return;
