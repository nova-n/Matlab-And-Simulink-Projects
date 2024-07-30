clc;
clear;

%A = [[1,9,8,1];[3,6,5,5];[7,6,5,8];[4,5,6,7]];
%B = [1;0;1;0];

%A = [[1,2];[2,4]];
%B = [1;2]

%A = [[0,4,-2];[6,-2,1];[4,8,-4]];
%B = [2;29;24];

A = [[-2,-17,4,3];[7,0,3,-2];[0,2,8,-6];[5,-13,-1,5]]
B = [0;0;-20;16]

A^-1

disp("Actual Answer:")
disp((A^-1)*B)

doneCondition = false;

A
%%Will first make sure top row has a leading number
if A(1,1) == 0
     holdingA = A(1,:);
     holdingB = B(1,:);
     nonLeadingZeroesRows = find(A(:,1) ~=0);
     %swap with the last row with a non leading 0
     A(1,:) = A(nonLeadingZeroesRows(end),:)
     B(1,:) = B(nonLeadingZeroesRows(end),:)
     A(nonLeadingZeroesRows(end),:) = holdingA;
     B(nonLeadingZeroesRows(end),:) = holdingB;
     disp( strcat("P(1," + num2str(nonLeadingZeroesRows(end)) + ")" ) )
     A
     B
end
disp( strcat("M_1(" + num2str(1/A(1,1)),")" ) )
B(1,:) = 1/A(1,1) * B(1,:);
A(1,:) = 1/A(1,1) * A(1,:)
B
%%Gets into triangular form
for i = 1:length(A) - 1 % i is the column of interest. Want to ignore the last column
    %will make all leading numbers except first row to be 0, and so on...
    for ii = 1:height(A)-1
        currentRow = height(A)- (ii - 1);% 0 to height(A)-1 would do the entire, from bottom to top, so subtract again
        if A(currentRow,i) ~=0 && currentRow ~= i
            %disp(A(currentRow,:))
            fac = -1*A(currentRow,i)/A(i,i);
            disp(strcat("A_", num2str(i) , " to " , num2str(currentRow) ," (" , num2str(fac) , ")" ) )
            B(currentRow,:) = fac * B(i,:) + B(currentRow,:);
            A(currentRow,:) = fac * A(i,:) + A(currentRow,:);
            A
            B
            arr = find(abs(A(currentRow,:)) > 10^-3); %cant use ~=0, since floating point.
            if length(arr) == 0
                continue;
            end
            mult = A(currentRow,arr(1) );
            disp(strcat("M_", num2str(currentRow) ," (1/" , num2str(mult) , ")" ) )
            B(currentRow,:) = B(currentRow,:) * 1/mult;
            A(currentRow,:) = A(currentRow,:) * 1/mult;
            A
            B
        end
    end
end
%%Will find any rows of all zeroes, and deletes them for now
z = [];
for i = 1:height(A)
    %need a tolerance since doing floating point
    if abs(sum(A(i,:)) ) <= 10^-3 
        z = [z,i]
    end
end
if length(z) ~=0
    for i = 1:length(z)
        if abs( sum(B(z(i) - (i-1) ,:) ) ) >=10^-3
            disp("system is inconsistent")
        end
        A(z(i) - (i-1) ,:) = [];
        B(z(i) - (i-1) ,:) = [];
    end
    % for i = 1:length(z)
    %     A(height(A) + i,:) = zeros(1,length(A));
    % end
    disp("removed zero row")
    A
    B
end

%%Gets into RREF
for i = 2:height(A)
    i;
   for ii = 1:i-1
       rowToAddTo =  ii;
       fac = -1*A(ii,i)/A(i,i);
       disp(strcat("A_", num2str(i) , " to " , num2str(ii) ," (" , num2str(fac) , ")" ) )
       A(rowToAddTo,:) = A(rowToAddTo,:) + fac*A(i,:);
       B(rowToAddTo,:) = B(rowToAddTo,:) + fac*B(i,:);
       A
   end
end
B