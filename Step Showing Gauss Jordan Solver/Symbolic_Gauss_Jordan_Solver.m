clc;
clear;

syms s

sI_A = sym(  [[s+2,-8,-7];[-5,s+4,-2];[9,3,s+1]]  );
Ident = sym(  [[1,0,0];[0,1,0];[0,0,1]]  );

%%testing
% sum( sI_A(1,:));
% symvar( sum( sI_A(1,:)) );
% symvar( sum( Ident(1,:)) );
% isempty( symvar( sum( Ident(1,:)) ) );

%%function call and work
%format rat
% [~,sI_A_Inv] = gaussJordanN(sI_A,Ident,"Problem 4");
% B = [1;5;-1];
% C = [7,1,2];
% transferFunction = simplify( simplifyFraction( C * sI_A_Inv * B ) );


%%function call and work
format rat
[~,sI_A_Inv] = gaussJordanN(sI_A,Ident,"Problem 3");
B = [0;1;0];
C = [1,0,1];
transferFunction = simplify( simplifyFraction( C * sI_A_Inv * B ) );

"\left[sI-A\right]^{-1} = " + latex(sI_A_Inv)
"$ \frac{Y\left(s\right)}{U\left(s\right)}=" + latex(sym(C))  + "" + latex(sI_A_Inv) + "" + latex(sym(B)) + "$"
"$ \frac{Y\left(s\right)}{U\left(s\right)}=" + latex(transferFunction) + "$"

% sI_A = sym(  [[s+1,-2,6];[4,s+5,0];[-3,3,s-7]]  );
% [~,sI_A_Inv] = gaussJordanN(sI_A,Ident,"Problem 4")
% B = [5;1;2];
% C = [2,0,2];
% transferFunction = simplifyFraction( C * sI_A_Inv * B )

% sI_A = sym(  [[0,-2,6];[4,s+5,0];[-3,3,s-7]]  );
% [~,sI_A_Inv] = gaussJordanN(sI_A,Ident,"Problem 4")
% B = [1;5;-1];
% C = [7,1,2];
% transferFunction = simplifyFraction( C * sI_A_Inv * B )

%disp("Actual Answer:");
%disp((A^-1)*B);

%%Will first make sure top row has a leading number
function [L,R] = gaussJordanN(L,R,problemName)
    disp("------   " + problemName + "   ------");
    L;
    R;
    "L = " + latex(L)
    "R = " + latex(R)
    if L(1,1) == 0
         holdingL = L(1,:);
         holdingR = R(1,:);
         nonLeadingZeroesRows = find(L(:,1) ~=0);
         %swap with the last row with a non leading 0
         L(1,:) = L(nonLeadingZeroesRows(end),:);
         R(1,:) = R(nonLeadingZeroesRows(end),:);
         L(nonLeadingZeroesRows(end),:) = holdingL;
         R(nonLeadingZeroesRows(end),:) = holdingR;
         disp( strcat("P(1," + num2str(nonLeadingZeroesRows(end)) + ") ... In Matlab: " , "L(1,:) = " , "L(" , num2str(nonLeadingZeroesRows(end)) , ",:) ; ",...
             "R(1,:) = " , "R(" , num2str(nonLeadingZeroesRows(end)) , ",:)") )
         
         L;
         R;
        "L = " + latex(L)
        "R = " + latex(R)
    end
    if L(1,1) ~=1
        disp( strcat("M_1(" + char(1/L(1,1)),") ... In Matlab: "  ,"R(1,:) = ","(1/L(1,1))", "*R(1,:)" ," ; L(1,:) = " ,"1/L(1,1))", "*L(1,:)") )
        R(1,:) = 1/L(1,1) * R(1,:);
        L(1,:) = 1/L(1,1) * L(1,:);
        "L = " + latex(L)
        "R = " + latex(R)
    end
    %%Gets into triangular form
    for i = 1:length(L) - 1 % i is the column of interest. Want to ignore the last column
        %will make all leading numbers except first row to be 0, and so on...
        if isempty( symvar( L(i,i) ) ) == 0 %eqiuvalent of checking if abs(fac) > 0      
            holdingL = L(i,:);
            holdingR = R(i,:);
            nonLeadingZeroesRows = find(L(:,i) ~=0);
            %swap with the last row with a non leading 0
            L(i,:) = L(nonLeadingZeroesRows(end),:);
            R(i,:) = R(nonLeadingZeroesRows(end),:);
            L(nonLeadingZeroesRows(end),:) = holdingL;
            R(nonLeadingZeroesRows(end),:) = holdingR;
            disp( strcat("P(" + num2str(i) + "," +  num2str(nonLeadingZeroesRows(end)) + ") ... In Matlab: " , "L(" , num2str(i) ,",:) = " , "L(" , num2str(nonLeadingZeroesRows(end)) , ",:) ; ",...
                "R(" , num2str(i) ,",:) = " , "R(" , num2str(nonLeadingZeroesRows(end)) , ",:)") )
            L = simplifyFraction(L);
            R = simplifyFraction(R);
            L;
            R;
            "L = " + latex(L)
            "R = " + latex(R)
        end
        for ii = 1:height(L)-1
            currentRow = height(L)- (ii - 1);% 0 to height(A)-1 would do the entire, from bottom to top, so subtract again
            if L(currentRow,i) ~=0
                if currentRow ~= i
                    %disp(L(currentRow,:))
                    fac = -1*L(currentRow,i)/L(i,i);
                    R(currentRow,:) = fac * R(i,:) + R(currentRow,:);
                    L(currentRow,:) = fac * L(i,:) + L(currentRow,:);
                    %if isempty( symvar(fac) ) == 0 %eqiuvalent of checking if abs(fac) > 0    
                        disp(strcat("A_", num2str(i) , " to " , "A_", num2str(currentRow) ," (" , char(fac) , ") ... In Matlab: " , "R(" , num2str(currentRow) , ",:) = " , "-1*L(",num2str(currentRow), ",", num2str(i),")"  , "* R(",num2str(i),",:)",...
                        " ; L(" , num2str(currentRow) , ",:) = " ,  "-1*L(",num2str(currentRow), ",", num2str(i),")" , "* L(",num2str(i),",:)") )
                        L = simplifyFraction(L);
                        R = simplifyFraction(R);
                        L;
                        R;
                        "L = " + latex(L)
                        "R = " + latex(R)
                    %end
                end
                arr = find(abs(L(currentRow,:)) > 10^-3); %cant use ~=0, since floating point.
                if length(arr) == 0
                    continue;
                end
                mult = L(currentRow,arr(1) );
                R(currentRow,:) = R(currentRow,:) * 1/mult;
                L(currentRow,:) = L(currentRow,:) * 1/mult;
                if mult ~=1
                    disp(strcat("M_", num2str(currentRow) ," (1/" , char(mult) , ") ... In Matlab: "  , "R(",num2str(currentRow),",:) = ", "(1/" , "L(",num2str(currentRow),",",char(arr(1)) , "))*R(" ,num2str(currentRow),",:)" ,...
                    " ; L(",num2str(currentRow),",:) = ", "(1/" , "L(",num2str(currentRow),",",char(arr(1)) , "))*L(" ,num2str(currentRow),",:)") )
                    L;
                    R;
                    "L = " + latex(L)
                    "R = " + latex(R)
                end
            end
        end
    end
    %%Will find any rows of all zeroes, and deletes them for now
    z = [];
    for i = 1:height(L)
        %need a tolerance since doing floating point
        if isempty( symvar( sum( L(i,:)) ) ) == 1 %checks if no symbols left in row
            if abs( double ( sum( L(i,:)) ) ) <= 10^-3 
                z = [z,i];
            end
        end
    end
    if length(z) ~=0
        for i = 1:length(z)
            if abs( sum(R(z(i) - (i-1) ,:) ) ) >=10^-3
                disp("system is inconsistent")
            end
            L(z(i) - (i-1) ,:) = [];
            R(z(i) - (i-1) ,:) = [];
        end
        % for i = 1:length(z)
        %     L(height(A) + i,:) = zeros(1,length(A));
        % end
        disp("removed zero row")
        L = simplifyFraction(L);
        R = simplifyFraction(R);
        L;
        R;
        "L = " + latex(L)
        "R = " + latex(R)
    end
    
    %%Gets into RREF
    for i = 2:height(L)
        i;
       for ii = 1:i-1
           rowToAddTo =  ii;
           fac = -1*L(ii,i)/L(i,i);
           L(rowToAddTo,:) = L(rowToAddTo,:) + fac*L(i,:);
           R(rowToAddTo,:) = R(rowToAddTo,:) + fac*R(i,:);
           if isempty( symvar(fac) ) == 0 %eqiuvalent of checking if abs(fac) > 0        
               disp(strcat("A_", num2str(i) , " to " ,"A_", num2str(ii) ," (" , char(fac) , ") ... In Matlab: " , "R(",num2str(ii),",:) = " ,"-1*L(",num2str(ii),",",num2str(i),")/L(",num2str(i),",",num2str(i),")" ,"*R(",num2str(ii),",:)",...
                   " ; L(",num2str(ii),",:) = " ,"-1*L(",num2str(ii),",",num2str(i),")/L(",num2str(i),",",num2str(i),")" ,"*L(",num2str(ii),",:)") )
             L = simplifyFraction(L);
             R = simplifyFraction(R);
             L;
             R;
             "L = " + latex(L)
             "R = " + latex(R)
           end
       end
    end
    L;
    R;
    disp(newline)
end