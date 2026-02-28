function append_pdfs_RZ(output, varargin)
% APPEND_PDFS_RZ  Concatenate multiple PDFs into one, silently.
% Includes a dummy last-page hack to avoid Ghostscript dropping
% the final input PDF in some installations.

    if nargin < 2
        error('append_pdfs_RZ:NotEnoughInputs', ...
              'Usage: append_pdfs_RZ(output, input1, input2, ...)');
    end

    % Determine append vs create
    if exist(output, 'file') == 2
        tmpOutput  = [tempname '.pdf'];
        inputFiles = [{output} varargin];
    else
        tmpOutput  = output;
        inputFiles = varargin;
    end

    % ----- Create dummy page to protect the last real file -----
    dummyPDF = [tempname '_dummy.pdf'];
    hf = figure('Visible','off');
    axis off
    set(hf,'PaperPosition',[0 0 8.5 11])
    print(hf, dummyPDF, '-dpdf', '-fillpage');
    close(hf);

    inputFiles = [inputFiles, {dummyPDF}];
    % -----------------------------------------------------------

    % Ghostscript command file
    cmdfile = [tempname '.txt'];
    fh = fopen(cmdfile,'w');

    % NOTE: -q makes Ghostscript completely silent
    fprintf(fh, ...
        '-q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile="%s" -f', ...
        tmpOutput);

    % Add all inputs
    for k = 1:numel(inputFiles)
        fprintf(fh,' "%s"', inputFiles{k});
    end
    fclose(fh);

    % ---- Silent Ghostscript call ----
    try
        status = evalc('ghostscript(''@"'' + string(cmdfile) + ''"'' )'); %#ok<NASGU>
    catch
        % ignore output completely
    end
    % ---------------------------------

    % Cleanup
    if exist(cmdfile,'file'), delete(cmdfile); end
    if exist(dummyPDF,'file'), delete(dummyPDF); end

    % Move temp into final if needed
    if ~strcmp(tmpOutput, output)
        movefile(tmpOutput, output, 'f');
    end
end
