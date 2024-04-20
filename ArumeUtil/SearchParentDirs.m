function [path,foundpath] = SearchParentDirs(currdir,strmatch,maxits)

    % searches all the parent directories above currdir for a match. If no
    % match, it returns 0
    foundpath = 0;

    for i = 1:maxits
        path = fileparts(currdir);
        dircontents = dir(path);
        dircontents = {dircontents.name};
        if any(strcmp(dircontents,strmatch))
            foundpath = 1;
            break
        else
            currdir = path;
        end
    end

end