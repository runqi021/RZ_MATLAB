function R = imagej_roi_read(roiPath, H, W)
%IMAGEJ_ROI_READ  Read an ImageJ/Fiji .roi file (or a RoiSet .zip) into MATLAB.
%
%   R = imagej_roi_read(roiPath)          parse only
%   R = imagej_roi_read(roiPath, H, W)    also build R(k).mask for an H x W image
%
%   roiPath : a single '*.roi', a 'RoiSet*.zip', or a folder (all *.roi in it)
%
%   R is a struct array with fields:
%     name      source file name
%     type      ImageJ type code
%     typeName  'polygon'|'rect'|'oval'|'line'|'freeline'|'polyline'|'freehand'|
%               'traced'|'angle'|'point'|'noRoi'
%     bounds    [x y w h], 1-based, x = column, y = row
%     x, y      vertex coordinates, 1-based image coords (empty for rect/oval)
%     mask      logical H x W (only when H,W are given)
%
%   ImageJ stores 0-based coordinates with x = column, y = row; everything
%   returned here is converted to MATLAB's 1-based convention.
%
%   Area types (rect, oval, polygon, freehand, traced) produce a filled mask.
%   Line types produce no area -- mask is left empty and a warning is issued,
%   because averaging pixels "inside" a line is not defined.
%
%   Runqi Zhang / 2026-08-12

    arguments
        roiPath (1,:) char
        H double = []
        W double = []
    end

    files = resolve_inputs(roiPath);
    assert(~isempty(files), 'No .roi files found at: %s', roiPath);

    R = struct('name',{},'type',{},'typeName',{},'bounds',{},'x',{},'y',{},'mask',{});
    for i = 1:numel(files)
        r = parse_one(files{i});
        if ~isempty(H) && ~isempty(W)
            r.mask = roi_to_mask(r, H, W);
        end
        R(end+1) = r; %#ok<AGROW>
    end
end

% =========================================================================
function files = resolve_inputs(p)
    files = {};
    if isfolder(p)
        d = dir(fullfile(p,'*.roi'));
        files = arrayfun(@(a) fullfile(a.folder,a.name), d, 'uni', 0);
        return;
    end
    [~,~,ext] = fileparts(p);
    if strcmpi(ext,'.zip')
        % RoiSet.zip -- unpack to a temp dir and take every .roi inside
        tmp = tempname; mkdir(tmp);
        unzip(p, tmp);
        d = dir(fullfile(tmp,'**','*.roi'));
        files = arrayfun(@(a) fullfile(a.folder,a.name), d, 'uni', 0);
    elseif isfile(p)
        files = {p};
    end
end

% =========================================================================
function r = parse_one(f)
    fid = fopen(f,'r','ieee-be');            % ImageJ .roi is BIG-endian
    assert(fid > 0, 'Cannot open %s', f);
    c = onCleanup(@() fclose(fid));

    magic = fread(fid, 4, '*char').';
    assert(strcmp(magic,'Iout'), '%s is not an ImageJ ROI (magic "%s")', f, magic);

    version = fread(fid, 1, 'int16');
    type    = fread(fid, 1, 'uint8');
    fread(fid, 1, 'uint8');                  % reserved

    top    = fread(fid, 1, 'int16');
    left   = fread(fid, 1, 'int16');
    bottom = fread(fid, 1, 'int16');
    right  = fread(fid, 1, 'int16');
    nCoord = fread(fid, 1, 'uint16');

    fseek(fid, 50, 'bof');
    options = fread(fid, 1, 'int16');

    names = {'polygon','rect','oval','line','freeline','polyline', ...
             'noRoi','freehand','traced','angle','point'};
    if type >= 0 && type <= 10
        typeName = names{type+1};
    else
        typeName = sprintf('unknown(%d)', type);
    end

    x = []; y = [];
    if nCoord > 0
        SUB_PIXEL = 128;
        subpix = version >= 222 && bitand(options, SUB_PIXEL) ~= 0;

        fseek(fid, 64, 'bof');
        xi = fread(fid, nCoord, 'int16');
        yi = fread(fid, nCoord, 'int16');

        if subpix
            % float arrays follow the int16 ones and are absolute, not relative
            xf = fread(fid, nCoord, 'float32');
            yf = fread(fid, nCoord, 'float32');
            if numel(xf) == nCoord && numel(yf) == nCoord && all(isfinite(xf))
                x = xf; y = yf;
            end
        end
        if isempty(x)
            x = double(xi) + left;           % int16 coords are relative to the bbox
            y = double(yi) + top;
        end
        x = x + 1;  y = y + 1;               % 0-based -> 1-based
    end

    [~, nm, ex] = fileparts(f);
    r = struct('name', [nm ex], 'type', type, 'typeName', typeName, ...
               'bounds', [left+1, top+1, right-left, bottom-top], ...
               'x', x(:), 'y', y(:), 'mask', []);
end

% =========================================================================
function m = roi_to_mask(r, H, W)
    m = false(H, W);
    bx = r.bounds(1); by = r.bounds(2); bw = r.bounds(3); bh = r.bounds(4);

    switch r.typeName
        case 'rect'
            cols = max(1,bx) : min(W, bx+bw-1);
            rows = max(1,by) : min(H, by+bh-1);
            m(rows, cols) = true;

        case 'oval'
            % ellipse inscribed in the bounding box, tested at pixel centres
            [cc, rr] = meshgrid(1:W, 1:H);
            cx = bx + bw/2 - 0.5;
            cy = by + bh/2 - 0.5;
            m = ((cc - cx)/(bw/2)).^2 + ((rr - cy)/(bh/2)).^2 <= 1;

        case {'polygon','freehand','traced'}
            assert(~isempty(r.x), 'ROI %s has no vertices', r.name);
            m = poly2mask(r.x, r.y, H, W);

        case 'point'
            idx = sub2ind([H W], round(r.y), round(r.x));
            m(idx) = true;

        otherwise
            warning('imagej_roi_read:noArea', ...
                ['ROI "%s" is type "%s", which encloses no area -- no mask built. ' ...
                 'Redraw it as a freehand/oval/polygon ROI if you want a dF/F trace.'], ...
                 r.name, r.typeName);
    end
end
