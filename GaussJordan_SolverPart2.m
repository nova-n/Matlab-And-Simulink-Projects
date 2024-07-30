clc;
clear;

A = [[1,9,8,1];[3,6,5,5];[7,6,5,8];[4,5,6,7]];
B = [1;0;1;0];

disp("Actual Answer:")
disp((A^-1)*B)

doneCondition = false;

A
%Will first make sure top row has a leading number
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
for i = 1:length(A) - 1 % i is the column of interest. Want to ignore the last column
    %will make all leading numbers except first row to be 0, and so on...
    for ii = 1:height(A)-1
        currentRow = height(A)- (ii - 1)% 0 to height(A)-1 would do the entire, from bottom to top, so subtract again
        if A(currentRow,i) ~=0 && currentRow ~= i
            %disp(A(currentRow,:))
            disp(strcat("A_", num2str(i) , " to " , num2str(currentRow) ," (" + num2str(-1*A(currentRow,i)/A(i,i)) + ")" ) )
            A(currentRow,1);
            A(currentRow,:) = -1*A(currentRow,i)/A(i,i) * A(i,:) +A(currentRow,:);
            A
        end
    end
end
%while doneConditoon == false
    
%end