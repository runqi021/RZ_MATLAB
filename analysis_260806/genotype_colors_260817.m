function C = genotype_colors_260817()
%GENOTYPE_COLORS_260817  The ONE genotype -> colour table.
%
%   C = genotype_colors_260817()   % struct, fieldname = genotype
%
% Every figure that colours by genotype reads this: the cartoon map's anatomy
% patches, the map's cell markers, and the group-activity overlays. It exists
% because those had drifted apart -- the map tinted Vglut2 anatomy [0.10 0.65
% 0.20] while a parallel table gave Vglut2 dots a lemon green, so one genotype
% appeared as two colours across a pair of figures meant to be read together.
%
% CHANGE A COLOUR HERE AND NOWHERE ELSE. A second hardcoded copy is how the
% drift happened the first time.
%
% Runqi Zhang / 2026-08-17

C = struct( ...
    'IO',     [0.50 0.50 0.50], ...   % grey
    'ChAT',   [0.85 0.10 0.10], ...   % red
    'Vglut2', [0.10 0.65 0.20], ...   % green
    'Vgat',   [0.45 0.65 0.95], ...   % light blue
    'Sst',    [0.55 0.20 0.75], ...   % purple
    'Sert',   [0.90 0.45 0.10]);      % orange
end
