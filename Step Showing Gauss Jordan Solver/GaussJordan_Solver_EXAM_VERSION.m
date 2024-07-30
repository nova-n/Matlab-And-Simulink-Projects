clc;
clear;

%A = [[1,9,8,1];[3,6,5,5];[7,6,5,8];[4,5,6,7]];
%B = [1;0;1;0];

%A = [[1,2];[2,4]];
%B = [1;2]

%A = [[0,4,-2];[6,-2,1];[4,8,-4]];
%B = [2;29;24];

%A = [[-2,-17,4,3];[7,0,3,-2];[0,2,8,-6];[5,-13,-1,5]]
%B = [0;0;-20;16]

%%Test Problems:

%2A
A = [[1,2,-1,1];[-1,1,2,-1];[2,-1,2,2];[1,1,-1,2]];
B = [6;3;14;8];
gaussJordanN(A,B,"Problem 2A")
%2B
A = [[1,2,-3,4];[2,2,-2,3];[0,1,1,0];[1,-1,1,-2]];
B = [12;10;-1;-4];
gaussJordanN(A,B,"Problem 2B")

%3A
A = [[2,0,1];[1,5,1];[-1,4,0]];
B = [[1,0,0];[0,1,0];[0,0,1]];
gaussJordanN(A,B,"Problem 3A")
%3B
A = [[29,-11,10];[-160,61,-55];[55,-21,19]];
B = [[1,0,0];[0,1,0];[0,0,1]];
gaussJordanN(A,B,"Problem 3B")

%4
A = [[0,0,-7,1];[0,0,5,0];[-7,5,0,2];[1,0,2,0]];
B = [0;0;0;0]; %for finding rowspace basis, no b vector needed. just set to 0 so not inconsistent system
gaussJordanN(A,B,"Problem 4")

%5
A = [[85,-28,-28];[10,-11,-11];[-46,-2,-2]];
B = [0;0;0];
gaussJordanN(A,B,"Problem 5 Eig Value = 0")

A = [[85 + 49 - sqrt(833)*i,-28,-28];[10,-11+ 49 - sqrt(833)*i,-11];[-46,-2,-2+ 49 - sqrt(833)*i]];
B = [0;0;0];
gaussJordanN(A,B, "Problem 5 Eig Value = 49 + sqrt(833)i")

A = [[85 - 49 - sqrt(833)*i,-28,-28];[10,-11- 49 - sqrt(833)*i,-11];[-46,-2,-2- 49 - sqrt(833)*i]];
B = [0;0;0];
gaussJordanN(A,B,"Problem 5 Eig Value = 49 - sqrt(833)i")

%disp("Actual Answer:");
%disp((A^-1)*B);

%%Will first make sure top row has a leading number
function gaussJordanN(A,B,problemName)
    disp("------   " + problemName + "   ------");
    A
    B
    if A(1,1) == 0
         holdingA = A(1,:);
         holdingB = B(1,:);
         nonLeadingZeroesRows = find(A(:,1) ~=0);
         %swap with the last row with a non leading 0
         A(1,:) = A(nonLeadingZeroesRows(end),:);
         B(1,:) = B(nonLeadingZeroesRows(end),:);
         A(nonLeadingZeroesRows(end),:) = holdingA;
         B(nonLeadingZeroesRows(end),:) = holdingB;
         disp( strcat("P(1," + num2str(nonLeadingZeroesRows(end)) + ") ... In Matlab: " , "A(1,:) = " , "A(" , num2str(nonLeadingZeroesRows(end)) , ",:) ; ",...
             "B(1,:) = " , "B(" , num2str(nonLeadingZeroesRows(end)) , ",:)") )
         A
         B
    end
    if A(1,1) ~=1
        disp( strcat("M_1(" + num2str(1/A(1,1)),") ... In Matlab: "  ,"B(1,:) = ","(1/A(1,1))", "*B(1,:)" ," ; A(1,:) = " ,"1/A(1,1))", "*A(1,:)") )
        B(1,:) = 1/A(1,1) * B(1,:);
        A(1,:) = 1/A(1,1) * A(1,:)
        B
    end
    %%Gets into triangular form
    for i = 1:length(A) - 1 % i is the column of interest. Want to ignore the last column
        %will make all leading numbers except first row to be 0, and so on...
        if abs(A(i,i)) <= 10^-3
            holdingA = A(i,:);
            holdingB = B(i,:);
            nonLeadingZeroesRows = find(A(:,i) ~=0);
            %swap with the last row with a non leading 0
            A(i,:) = A(nonLeadingZeroesRows(end),:);
            B(i,:) = B(nonLeadingZeroesRows(end),:);
            A(nonLeadingZeroesRows(end),:) = holdingA;
            B(nonLeadingZeroesRows(end),:) = holdingB;
            disp( strcat("P(" + num2str(i) + "," +  num2str(nonLeadingZeroesRows(end)) + ") ... In Matlab: " , "A(" , num2str(i) ,",:) = " , "A(" , num2str(nonLeadingZeroesRows(end)) , ",:) ; ",...
                "B(" , num2str(i) ,",:) = " , "B(" , num2str(nonLeadingZeroesRows(end)) , ",:)") )
            A
            B
        end
        for ii = 1:height(A)-1
            currentRow = height(A)- (ii - 1);% 0 to height(A)-1 would do the entire, from bottom to top, so subtract again
            if A(currentRow,i) ~=0
                if currentRow ~= i
                    %disp(A(currentRow,:))
                    fac = -1*A(currentRow,i)/A(i,i);
                    B(currentRow,:) = fac * B(i,:) + B(currentRow,:);
                    A(currentRow,:) = fac * A(i,:) + A(currentRow,:);
                    if abs(fac)>0
                        disp(strcat("A_", num2str(i) , " to " , num2str(currentRow) ," (" , num2str(fac) , ") ... In Matlab: " , "B(" , num2str(currentRow) , ",:) = " , "-1*A(",num2str(currentRow), ",", num2str(i),")"  , "* B(",num2str(i),",:)",...
                        " ; A(" , num2str(currentRow) , ",:) = " ,  "-1*A(",num2str(currentRow), ",", num2str(i),")" , "* A(",num2str(i),",:)") )
                        A
                        B
                    end
                end
                arr = find(abs(A(currentRow,:)) > 10^-3); %cant use ~=0, since floating point.
                if length(arr) == 0
                    continue;
                end
                mult = A(currentRow,arr(1) );
                B(currentRow,:) = B(currentRow,:) * 1/mult;
                A(currentRow,:) = A(currentRow,:) * 1/mult;
                if mult ~=1
                    disp(strcat("M_", num2str(currentRow) ," (1/" , num2str(mult) , ") ... In Matlab: "  , "B(",num2str(currentRow),",:) = ", "(1/" , "A(",num2str(currentRow),",",num2str(arr(1)) , "))*B(" ,num2str(currentRow),",:)" ,...
                    " ; A(",num2str(currentRow),",:) = ", "(1/" , "A(",num2str(currentRow),",",num2str(arr(1)) , "))*A(" ,num2str(currentRow),",:)") )
                    A
                    B
                end
            end
        end
    end
    %%Will find any rows of all zeroes, and deletes them for now
    z = [];
    for i = 1:height(A)
        %need a tolerance since doing floating point
        if abs(sum(A(i,:)) ) <= 10^-3 
            z = [z,i];
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
           A(rowToAddTo,:) = A(rowToAddTo,:) + fac*A(i,:);
           B(rowToAddTo,:) = B(rowToAddTo,:) + fac*B(i,:);
           if abs(fac)>0         
               disp(strcat("A_", num2str(i) , " to " , num2str(ii) ," (" , num2str(fac) , ") ... In Matlab: " , "B(",num2str(ii),",:) = " ,"-1*A(",num2str(ii),",",num2str(i),")/A(",num2str(i),",",num2str(i),")" ,"*B(",num2str(ii),",:)",...
                   " ; A(",num2str(ii),",:) = " ,"-1*A(",num2str(ii),",",num2str(i),")/A(",num2str(i),",",num2str(i),")" ,"*A(",num2str(ii),",:)") )
             A
             B
           end
       end
    end
    A;
    B;
    disp(newline)
end