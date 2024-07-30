clc;
clear;

A = [[1,9,8];[3,6,5];[7,6,5]];
B = [1;0;1];

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
%will make all leading numbers except first row to be 0
for i = 0:height(A)-2 %-1 would do the entire, form bottom to top. Want to stop at first row.
    currentRow = height(A)-i;
    if A(currentRow,1) ~= 0
        %disp(A(currentRow,:))
        disp(strcat("A_", num2str(1) , " to " , num2str(currentRow) ," (" + num2str(-1*A(currentRow,1)/A(1,1)) + ")" ) )
        A(currentRow,1);
        A(currentRow,:) = -1*A(currentRow,1)/A(1,1) * A(1,:) +A(currentRow,:);
        A
    end
end

%while doneConditoon == false
    
%end