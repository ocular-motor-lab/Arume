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


% Now do the same for the two additional paths we need: FilterConfigFiles
% and Shaders
rootdir = fileparts(new_arume_folder);
filelist = dir(fullfile(rootdir, ['**',filesep,'*']));  % get list of files and folders in any subfolder
filelist = filelist([filelist.isdir]);  % remove folders from list

% add shaders dir
shadersdirs = filelist(contains({filelist.name},'Shaders'));

% remove all existing paths, ignore warnings:
warning('off')
for i = 1:length(shadersdirs)
    fname = [shadersdirs(i).folder,filesep,shadersdirs(i).name];
    disp(['REMOVED: ',fname])
    rmpath(fname)

    fname = [shadersdirs(i).folder,filesep,'FilterConfigFiles'];
    disp(['REMOVED: ',fname])
    rmpath(fname)

    fname = [shadersdirs(i).folder,filesep,'VisualStimuli'];
    disp(['REMOVED: ',fname])
    rmpath(fname)
end
warning('on');

% if nocla exists, use nocla. If not, choose the latest deployment version
isnocla = contains({shadersdirs.folder},'nocla');
if any(isnocla)
    if sum(isnocla) == 1
        shadersdir = shadersdirs(isnocla);
        addpath(shadersdir.folder);
        disp(['ADDED TO THE PATH ' shadersdir.folder])
    else
        error('Multiple nocla folders detected. Aborting setup.')
    end
else
    shadersdir = shadersdirs(end);
end

addpath([shadersdir.folder, filesep, shadersdir.name]);
disp(['ADDED TO THE PATH ' [shadersdir.folder, filesep, shadersdir.name]])

% and add filterconfig, which is always in the same parent folder
addpath([shadersdir.folder, filesep, 'FilterConfigFiles']);
disp(['ADDED TO THE PATH ' [shadersdir.folder, filesep, 'FilterConfigFiles']])

addpath([shadersdir.folder, filesep, 'VisualStimuli']);
disp(['ADDED TO THE PATH ' [shadersdir.folder, filesep, 'VisualStimuli']])

savepath;
