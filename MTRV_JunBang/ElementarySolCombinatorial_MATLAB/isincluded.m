function [inclusion,ind] = isincluded(col,col_mat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check the Inclusion of the Column                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% col: column (debris sequence)
% col_mat: column matrix
% inclusion = 0: col is not included in the matrix
%             1: col is aleady included in the matrix
% ind: index of row in the matrix

if(~isempty(col_mat))
    inclusion = 0;
    if(iscell(col_mat))
        for ind=1:length(col_mat(:,1))
            if(length(col_mat{ind,1})==length(col) && sum(col_mat{ind,1}==col)==length(col))
                inclusion = 1;
                break;
            end
        end
    else
        for ind=1:length(col_mat)
            if(col==col_mat(ind))
                inclusion = 1;
                break;
            end
        end
    end
else
    inclusion = 0;
    ind = 0;
end
end