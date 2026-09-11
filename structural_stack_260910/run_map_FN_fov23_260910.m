% run_map_FN_fov23_260910.m
% -----------------------------------------------------------------------
%  Prep + render map_FN_fov2 and fov3, the other two tiles of the 260909
%  facial-nucleus strip.  Same parameters as fov1 -- the three tiles differ
%  only in stage x (-1727.12 / -1330.11 / -933.11 um, step 397.0), so anything
%  that differed between them would be a bug, not a choice.
%
%  Runqi Zhang / 2026-09-10
% -----------------------------------------------------------------------

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)); addpath(here);

ROOT     = 'C:\Users\Admin\Desktop\260909_ChAT_g8m_Shiverer';
STEMS    = {'map_FN_fov2_00001', 'map_FN_fov3_00001'};
ZSURFACE = 0;

for i = 1:numel(STEMS)
    stem = STEMS{i};
    tif  = fullfile(ROOT, stem, [stem '.tif']);
    fprintf('\n\n############ %s ############\n', stem);
    struct_stack_prep_260910(tif, 'PixTarget', 1.0, 'ZSurface', ZSURFACE);
    struct_stack_render_260910(fullfile(ROOT, stem, [stem '_stackmeta.mat']), ...
                               'BinUm', 100, 'ClimPct', [0.5 99.5]);
end
