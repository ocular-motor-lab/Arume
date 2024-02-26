function dataTable = VertCatTablesMissing(t1,t2)
% Concatenates two tables that may have different variables
% adding missing elements where it corresponds.
%
% Example: 
%
% t1 =
% 
%   3×3 table
% 
%     Var1      Var2      Var6 
%     ____    ________    _____
% 
%       1     "one"       true 
%      20     "twenty"    false
%     300     "three"     false
% 
% t2
% 
% t2 =
% 
%   3×5 table
% 
%       Var3      Var4     Var5      Var1    Var7
%     ________    ____    _______    ____    ____
% 
%     {[2 43]}     as     someone      3      3  
%     {[7 73]}     ba     then         5      4  
%     {[6 63]}     ca     test        66      5  
% 
% VertCatTablesMissing(t1,t2)
% 
% ans =
% 
%   6×7 table
% 
%         Var3        Var4    Var1      Var2       Var6       Var5        Var7
%     ____________    ____    ____    _________    ____    ___________    ____
% 
%     {0×0 double}              1     "one"          1     <undefined>    NaN 
%     {0×0 double}             20     "twenty"       0     <undefined>    NaN 
%     {0×0 double}            300     "three"        0     <undefined>    NaN 
%     {[    2 43]}     as       3     <missing>    NaN     someone          3 
%     {[    7 73]}     ba       5     <missing>    NaN     then             4 
%     {[    6 63]}     ca      66     <missing>    NaN     test             5 
   
% 
%     
if ( ~istable(t1) )
    t1 = table();
end
if ( ~istable(t2) )
    t2 = table();
end

% find logical columns and make them doubles so we can add nans in the
% missing rows

for colnamecell = t1.Properties.VariableNames
    colname = colnamecell{1};
    if ( islogical(t1.(colname) ) )
        t1.(colname) = double(t1.(colname));
    end
end
for colnamecell = t2.Properties.VariableNames
    colname = colnamecell{1};
    if ( islogical(t2.(colname)) )
        t2.(colname) = double(t2.(colname));
    end
end


% find the columns in table 2 that are missing in table 1
t1colmissing = setdiff(t2.Properties.VariableNames, t1.Properties.VariableNames,'stable');
% find the columns in table 1 that are missing in table 2
t2colmissing = setdiff(t1.Properties.VariableNames, t2.Properties.VariableNames,'stable');

% cannot add missing to cell columns or char columns 
% so here we find the missing columns that are cells and chars
% first in table 1
t1colmissingCells = {};
t1colmissingChars = {};
for colname = t1colmissing
    if iscell(t2.(colname{1}))
        t1colmissingCells(end+1) = colname;
    end
    if ischar(t2.(colname{1}))
        t1colmissingChars(end+1) = colname;
    end
end

% Then find the cell and char missing columns in in table 2
t2colmissingCells = {};
t2colmissingChars = {};
for colname = t2colmissing
    if iscell(t1.(colname{1}))
        t2colmissingCells(end+1) = colname;
    end
    if ischar(t1.(colname{1}))
        t2colmissingChars(end+1) = colname;
    end
end

% dealing with cell columns by adding empty cells
for colname = t1colmissingCells
    t1.(colname{1}) = cell(height(t1), 1);
end
for colname = t2colmissingCells
    t2.(colname{1}) = cell(height(t2), 1);
end

% dealing with char columns add empty spaces. Does not work with empty
% strings
for colname = t1colmissingChars
    t1.(colname{1}) = repmat(' ',height(t1), width(t2.(colname{1})));
end
for colname = t2colmissingChars
    t2.(colname{1}) = repmat(' ',height(t2), width(t1.(colname{1})));
end

% We've added the columns so they are not missing anymore
t1colmissing = setdiff(t1colmissing, t1colmissingCells,'stable');
t1colmissing = setdiff(t1colmissing, t1colmissingChars,'stable');
t2colmissing = setdiff(t2colmissing, t2colmissingCells,'stable');
t2colmissing = setdiff(t2colmissing, t2colmissingChars,'stable');

% find the common columns that are not missing in either table
t1andt2 = intersect(t2.Properties.VariableNames, t1.Properties.VariableNames, 'stable');

% first concatenate the common columns
t = [t1(:,t1andt2);t2(:,t1andt2)];
% then add to the unique columns of t1 missing
% rows for each element of t2
t12 = t1(:,t2colmissing);
if ( height(t2)>0 && ~isempty(t2colmissing))
    t12{height(t1)+(1:height(t2)),:} = missing;
end
% then add to the unique columns of t2 missing
% rows for each element of t1
t21 = t2(:,t1colmissing);
if ( height(t1)>0 && ~isempty(t1colmissing))
    t21{height(t2)+(1:height(t1)),:} = missing;
    % move the values up so the rows match
    t21 = [t21(height(t2)+1:end,:);t21(1:height(t2),:)];
end

% concatenate the 3 tables
if ( ~isempty( t) && ~isempty( t12) && ~isempty( t21) )
    dataTable = [t t12 t21];
elseif ( ~isempty( t) && ~isempty( t12) && isempty( t21) )
    dataTable = [t t12];
elseif ( ~isempty( t) && isempty( t12) && ~isempty( t21) )
    dataTable = [t t21];
elseif ( isempty( t) && ~isempty( t12) && ~isempty( t21) )
    dataTable = [t12 t21];
elseif ( ~isempty( t) && isempty( t12) && isempty( t21) )
    dataTable = t;
elseif ( isempty( t) && ~isempty( t12) && isempty( t21) )
    dataTable = t12;
elseif ( isempty( t) && isempty( t12) && ~isempty( t21) )
    dataTable = t21;
end

end