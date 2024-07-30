clc;
clear;

syms s

sI_A = sym(  [[s+2,-8,-7];[-5,s+4,-2];[9,3,s+1]]  );
Ident = sym(  [[1,0,0];[0,1,0];[0,0,1]]  );

%%testing
sum( sI_A(1,:))
symvar( sum( sI_A(1,:)) )
symvar( sum( Ident(1,:)) )
isempty( symvar( sum( Ident(1,:)) ) )

%%function call and work
format rat
[~,sI_A_Inv] = gaussJordanN(sI_A,Ident,"Problem 4")
B = [1;5;-1];
C = [7,1,2];
C * sI_A_Inv * B



%disp("Actual Answer:");
%disp((A^-1)*B);

%%Will first make sure top row has a leading number
function [A,B] = gaussJordanN(A,B,problemName)
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
         disp( strcat("P(1," + char(nonLeadingZeroesRows(end)) + ") ... In Matlab: " , "A(1,:) = " , "A(" , char(nonLeadingZeroesRows(end)) , ",:) ; ",...
             "B(1,:) = " , "B(" , char(nonLeadingZeroesRows(end)) , ",:)") )
         A
         B
    end
    if A(1,1) ~=1
        disp( strcat("M_1(" + char(1/A(1,1)),") ... In Matlab: "  ,"B(1,:) = ","(1/A(1,1))", "*B(1,:)" ," ; A(1,:) = " ,"1/A(1,1))", "*A(1,:)") )
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
            disp( strcat("P(" + char(i) + "," +  char(nonLeadingZeroesRows(end)) + ") ... In Matlab: " , "A(" , char(i) ,",:) = " , "A(" , char(nonLeadingZeroesRows(end)) , ",:) ; ",...
                "B(" , char(i) ,",:) = " , "B(" , char(nonLeadingZeroesRows(end)) , ",:)") )
            A = simplifyFraction(A);
            B = simplifyFraction(B);
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
                        disp(strcat("A_", char(i) , " to " , char(currentRow) ," (" , char(fac) , ") ... In Matlab: " , "B(" , char(currentRow) , ",:) = " , "-1*A(",char(currentRow), ",", char(i),")"  , "* B(",char(i),",:)",...
                        " ; A(" , char(currentRow) , ",:) = " ,  "-1*A(",char(currentRow), ",", char(i),")" , "* A(",char(i),",:)") )
                        A = simplifyFraction(A);
                        B = simplifyFraction(B);
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
                    disp(strcat("M_", char(currentRow) ," (1/" , char(mult) , ") ... In Matlab: "  , "B(",char(currentRow),",:) = ", "(1/" , "A(",char(currentRow),",",char(arr(1)) , "))*B(" ,char(currentRow),",:)" ,...
                    " ; A(",char(currentRow),",:) = ", "(1/" , "A(",char(currentRow),",",char(arr(1)) , "))*A(" ,char(currentRow),",:)") )
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
        if isempty( symvar( sum( A(i,:)) ) ) == 1 %checks if no symbols left in row
            if abs( double ( sum( A(i,:)) ) ) <= 10^-3 
                z = [z,i];
            end
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
        A = simplifyFraction(A);
        B = simplifyFraction(B);
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
           if isempty( symvar(fac) ) == 0 %eqiuvalent of checking if abs(fac) > 0        
               disp(strcat("A_", char(i) , " to " , char(ii) ," (" , char(fac) , ") ... In Matlab: " , "B(",char(ii),",:) = " ,"-1*A(",char(ii),",",char(i),")/A(",char(i),",",char(i),")" ,"*B(",char(ii),",:)",...
                   " ; A(",char(ii),",:) = " ,"-1*A(",char(ii),",",char(i),")/A(",char(i),",",char(i),")" ,"*A(",char(ii),",:)") )
             A = simplifyFraction(A);
             B = simplifyFraction(B);
             A
             B
           end
       end
    end
    A;
    B;
    disp(newline)
end