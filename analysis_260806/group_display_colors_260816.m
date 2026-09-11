function DISP = group_display_colors_260816()
%GROUP_DISPLAY_COLORS_260816  The one definition of the display groups + colours.
%
%   DISP = group_display_colors_260816()
%
% Nx2 cell: {display name, RGB}. Display name is "<genotype> <class>", exactly
% as built by  string(group) + " " + string(class)  from cell_classes_260816 --
% except IO, which is its own display group regardless of class.
%
% WHY THIS IS A FUNCTION and not a literal in each script: the population-average
% overlays and the cartoon map have to use the SAME colour for the same group, or
% the reader pairs a green trace with a green dot that is a different population.
% These colours have already been changed once by hand (2026-08-16: Vgat I and
% Vgat tonic swapped, Vglut2 pre-I moved to yellow-green and Vglut2 I to the dark
% green), and doing that in two places is how they drift apart.
%
% Runqi Zhang / 2026-08-16

DISP = { 'IO',                            [0    0    0   ]    % black
         'ChAT post-I',                   [0.85 0.10 0.10]    % red
         'Vglut2 pre-I',                  [0.62 0.76 0.10]    % yellow green
         'Vglut2 I',                      [0.05 0.45 0.15]    % dark green
         'Vgat I',                        [0.45 0.65 0.95]    % light blue
         'Vgat tonic pre-I suppressed',   [0.05 0.20 0.70]    % dark blue
         'Sst null',                      [0.55 0.20 0.75]    % purple
         'Sert post-I',                   [0.90 0.45 0.10] }; % orange
end
