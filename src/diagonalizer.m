function [ G, U, S ] = diagonalizer( A , m , n)
%DIAGNOLIZER Summary of this function goes here
%   Detailed explanation goes here
%   Assumption: column 1:m contain only integers
%               colunm m+1:n is the result of M(integer)*v(symbols)
%   

G = eye(n);
U = eye(n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diagonalize the submatrix  A(:,1:m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('======== making matrix diagonal ========== ');

for i=1:m

    str = sprintf('=== diagonal element: %d : ', i);
    disp(str);
    
    while(1)
        skip = 0;
        %% STEP 1: interchange for the smallest diagonal element at i
        % check row and interchange rows
        col = A(:,i);
        col(col==0) = inf;
        [p, idx] = min(abs(col));
        T = rowcol_itrchg(i, idx, n);
        A = T*A;
        G = T*G;
        % check column and interchange colums
        row = A(i,1:m);
        row(row==0) = inf;
        [p, idx] = min(abs(row));
        T = rowcol_itrchg(i, idx, n);
        A = A*T;
        U = U*T;
        
%         seq = seq*T;

        %% STEP 2: check divisibility on rows
        for j=i+1:n
           if( mod(A(j,i), A(i,i)) ~= 0)
                % A(j,i) = q*A(i,i) + r
                q = floor(A(j,i)/A(i,i));
                A(j,:) = A(j,:) - q*A(i,:);
                T = rowcol_add(q, 1, i, j, n);
                G = T*G;
                skip = 1;
                break;
           end
        end
        if(skip == 1) 
            continue; 
        end
        
        %% STEP 3: check divisibility on columns
         for j=i+1:m
           if( mod(A(i,j), A(i,i)) ~= 0)
                % A(i,j) = q*A(i,i) + r
                q = floor(A(i,j)/A(i,i));
                A(:,j) = A(:,j) - q*A(:,i);
                T = rowcol_add(q, 2, i, j, n);
                U = U*T;
                skip = 1;
                break;
           end
        end
        if(skip == 1) 
            continue; 
        end       
            
        %% STEP 4 : elimination
        % eliminate rows
        for j=i+1:n
            if( A(j,i) ~= 0)
                % A(j,i) = q*A(i,i)
                q = A(j,i)/A(i,i);
                A(j,:) = A(j,:) - q*A(i,:);
                T = rowcol_add(q, 1, i, j, n);
                G = T*G;     
            end
        end
        % eliminate cols
        for j=i+1:m
            if( A(i,j) ~= 0)
                % A(i,j) = q*A(i,i)
                q = A(i,j)/A(i,i);
                A(:,j) = A(:,j) - q*A(:,i);
                T = rowcol_add(q, 2, i, j, n);
                U = U*T;     
            end
        end
        
        %% finish current diagonal element
        break;
    end
    
end

S = A;

end

