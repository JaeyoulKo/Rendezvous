function [col_next,ind] = getnextcol(col,n,iscut)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get Next Column in Lexicographic Order                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% col: current column (debris sequence)                   ex) [1,2,3]
% col_next: next column
% n: number of debris
% iscut = 0: add additional debris in the current one     ex) [1,2,3,4]
%         1: change the final debris in the current one   ex) [1,2,4]

if(nargin<3)
    iscut = 0;
end

if(col==n)
    col_next = [];
else
    l = length(col);
    node = col(l);
    if(node<n)
        if(~iscut)
            col_next = [col,node+1];
        else
            col_next = [col(1:l-1),node+1];
        end
    else
        if(l==2)
            col_next = col(l-1)+1;
        else
            col_next = [col(1:l-2),col(l-1)+1];
        end
    end
end
ind = length(col_next)-length(col);
end