% set the matlab working directly to this folder and run this script from
% the command line
%
%   arume_setup
%

% remove current arume from path
new_arume_folder = fileparts(mfilename('fullpath'));

paths=regexpi(path,['[^;]*arume[^;]*;'],'match');
if ( length( paths) == 0 )
    paths=regexp(path,['[^;]*arume[^;]*;'],'match');
end
if ( length( paths) > 0 )
    
    disp('This folders will be removed from the path:');
    disp(paths')
    response = input('do you want to continue? (y/n)','s');
    if ( lower(response) ~= 'y')
        return
    end
    for p=paths
        s=char(p);
        rmpath(s);
    end
end

response = input('The new folders will be added to the path, continue? (y/n)','s');
if ( lower(response) ~= 'y')
    return
end
addpath(new_arume_folder, fullfile(new_arume_folder, 'ArumeUtil'));
disp(['ADDED TO THE PATH ' new_arume_folder ])
disp(['ADDED TO THE PATH ' fullfile(new_arume_folder, 'ArumeUtil') ])



savepath;
