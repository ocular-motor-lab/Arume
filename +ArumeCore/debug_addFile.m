%%
LinkedFiles = struct();


LinkedFiles = addExistingFile(LinkedFiles, "C:\Users\omlab\Documents\Arume\+ArumeCore", 'EyelinkFile', ...
    'C:\secure\ClaraDataRaw\MotionEllipsesTesting\A1MotionEllipses__000__Zdfff\A1MotionEllipses__000__Zdfff_20260422_133140.edf');
LinkedFiles = addExistingFile(LinkedFiles, "C:\Users\omlab\Documents\Arume\+ArumeCore", 'EyelinkFile', ...
    'C:\secure\ClaraDataRaw\MotionEllipsesTesting\A1MotionEllipses__000__Zdfff\A1MotionEllipses__000__Zdfff_20260422_133336.edf');



function LinkedFiles = addExistingFile(LinkedFiles, dataPath, fileTag, filePath)


[~,fileName, ext] = fileparts(filePath);

if ( ~isfield(LinkedFiles, fileTag) )
    LinkedFiles.(fileTag) = [fileName ext];
else
    if ~iscell(LinkedFiles.(fileTag))
        LinkedFiles.(fileTag) = {LinkedFiles.(fileTag)};
    end
    LinkedFiles.(fileTag) = vertcat( LinkedFiles.(fileTag), [fileName ext] );
end

end

