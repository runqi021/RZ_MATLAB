% run_map_FN_fov1_260910.m
% -----------------------------------------------------------------------
%  Driver: structural 3-D stack -> processed volume -> projection + 3-D figures
%  for 260909_ChAT_g8m_Shiverer / map_FN_fov1_00001.
%
%  Stage 1 reads the raw 16.8 GB stack ONCE (~5 min) and caches the 1 um
%  isotropic volumes; stage 2 only reads those, so re-plotting is cheap. Delete
%  the *_avg1um.tif / *_stackmeta.mat to force a re-read.
%
%  Runqi Zhang / 2026-09-10
% -----------------------------------------------------------------------

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)); addpath(here);

TIF = 'C:\Users\Admin\Desktop\260909_ChAT_g8m_Shiverer\map_FN_fov1_00001\map_FN_fov1_00001.tif';

% ZSurface is the MOTOR z of the pia, and depth = ZSurface - z.  This stack was
% started at motor +30 and run to -500, which is the usual "start above the
% surface" pattern, so 0 is the operator's own z zero.  It is NOT measured here:
% stage 1 prints the (top-bot) half-rise it finds so the assumption can be
% checked, and changes nothing on the strength of it.
ZSURFACE = 0;

meta = struct_stack_prep_260910(TIF, 'PixTarget', 1.0, 'ZSurface', ZSURFACE);

metaMat = fullfile(fileparts(TIF), 'map_FN_fov1_00001_stackmeta.mat');
struct_stack_render_260910(metaMat, 'BinUm', 100, 'ClimPct', [0.5 99.5]);
