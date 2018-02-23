% Stokes FMMs in R^3.
%
% Particle FMM routines.
%   stfmm3dpart - Stokes particle FMM in R^3.
%   sthfmm3dpart - Stokes half-space particle FMM in R^3.
%   lfmm3dpartquad - Laplace particle FMM in R^3.
%
% Direct evaluation routines.
%   st3dpartdirect - Stokes particle interactions in R^3.
%   sth3dpartdirect - Stokes half-space particle interactions in R^3.
%   l3dpartquaddirect - Laplace particle interactions in R^3.
%
% Field postprocessing routines.
%   st3stress - Evaluate stress tensor from p and \grad u.
%   st3strain - Evaluate strain tensor from \grad u.
%   st3traction - Evaluate traction vector from p, \grad u, and n.
%   st3stress2traction - Evaluate traction vector from stress tensor and n.
%
% Internal utility functions.
%   stfmm3dprini   - Initialize internal printing routines.
%
% Demos.
%
%   test_st3dpart_direct - test Stokes particle FMM vs direct routines.
%   test_sth3dpart_direct - test Stokes half-space particle FMM.
%   test_sth3dpart_direct2 - test Fortran MEX vs Matlab routines.
%                            Half space Stokes particle interactions. 
%   test_sth3dpart_direct3 - test Fortran MEX vs Matlab routines. 
%                            Half space Stokes particle FMM.
%


%% Copyright (C) 2009-2012: Leslie Greengard and Zydrunas Gimbutas
%% Contact: greengard@cims.nyu.edu
%% 
%% This program is free software; you can redistribute it and/or modify 
%% it under the terms of the GNU General Public License as published by 
%% the Free Software Foundation; either version 2 of the License, or 
%% (at your option) any later version.  This program is distributed in 
%% the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
%% even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
%% PARTICULAR PURPOSE.  See the GNU General Public License for more 
%% details. You should have received a copy of the GNU General Public 
%% License along with this program; 
%% if not, see <http://www.gnu.org/licenses/>.
