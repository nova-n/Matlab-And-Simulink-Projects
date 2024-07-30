clc;
clear;

A = [[1,9,8,1];[3,6,5,5];[7,6,5,8];[4,5,6,7]];
B = [1;0;1;0];

disp("Actual Answer:")
disp((A^-1)*B)

doneCondition = false;

A
%%Will first make sure top row has a leading number
if A(1,1) == 0
     holding = A(1,:);
     nonLeadingZeroesRows = find(A(:,1) ~=0);
     %swap with the last row with a non leading 0
     A(1,:) = A(nonLeadingZeroesRows(end),:)
     A(nonLeadingZeroesRows(end),:) = holding;
     disp( strcat("P(1," + num2str(nonLeadingZeroesRows(end)) + ")" ) )
     A
end
disp( strcat("M_1(" + num2str(1/A(1,1)),")" ) )
A(1,:) = 1/A(1,1) * A(1,:)
%%Gets into triangular form
for i = 1:length(A) - 1 % i is the column of interest. Want to ignore the last column
    %will make all leading numbers except first row to be 0, and so on...
    for ii = 1:height(A)-1
        currentRow = height(A)- (ii - 1)% 0 to height(A)-1 would do the entire, from bottom to top, so subtract again
        if A(currentRow,i) ~=0 && currentRow ~= i
            %disp(A(currentRow,:))
            fac = -1*A(currentRow,i)/A(i,i);
            disp(strcat("A_", num2str(i) , " to " , num2str(currentRow) ," (" , num2str(fac) , ")" ) )
            A(currentRow,:) = fac * A(i,:) + A(currentRow,:);
            A
            arr = find(A(currentRow,:) ~= 0);
            mult = A(currentRow,arr(1) );
            disp(strcat("M_", num2str(currentRow) ," (1/" , num2str(mult) , ")" ) )
            A(currentRow,:) = A(currentRow,:) * 1/mult;
            A
        end
    end
end
%%Will find any rows of all zeroes, and move them to the bottom
z = [];
for i = 1:height(A)
    %need a tolerance since doing floating point
    z = [z,find( abs(sum(A(i,:)) - 1e-3 ) == 0 )];
end
if length(z) ~=0
    for i = 1:length(z)
        A(z(i) - (i-1),:) = [];
    end
    for i = 1:length(z)
        A(height(A) + i,:) = zeros(1,length(A));
    end
end
%%

for i = 2:height(A)
    i;
   for ii = 1:i-1
       rowToAddTo =  ii;
       fac = -1*A(ii,i)/A(i,i);
       disp(strcat("A_", num2str(i) , " to " , num2str(ii) ," (" , num2str(fac) , ")" ) )
       A(rowToAddTo,:) = A(rowToAddTo,:) + fac*A(i,:);
       A
   end
end