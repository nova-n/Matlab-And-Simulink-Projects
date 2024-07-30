clc;
clear;
close all;

digitsDesired = 3; %%%enter digits here!
primes = primesUpToDigit(digitsDesired) %%%call it here!!!
sprintf("Sum of primes up to %d digits: %d",digitsDesired, sum(primes))

function [primes] = primesUpToDigit(digits)
    maxNum = 10^digits - 1; %for example, if want 3 digits, then get 10^3 -1 = 1000 -1 = 999, the largest 3 digit numbner
    primes = [2,3]; %will get populated as the loop progresses
    startNumber = primes(end)+1; %will start collecting at the end of the primes list + 1, so in this case, at 3+1=4
    
    for i = startNumber:maxNum
        for ii = 1:length(primes)
            if mod(i,primes(ii)) == 0 %modulo to check for remainders of dividing i by primes(ii)
                %if remainder is 0, then that number is divisible, so break
                %loop, and skip this number
                break
            end
        end
        %if it reached the end of the list without breaking loop (for ii), then
        %that number was not found to be divisible by anything, so its prime
        %this works because you keep adding more and more numbers to the list,
        %so bigger numbers have more factors to check against
        if ii == length(primes) 
            primes(end+1) = i;
        end
    end
    primes;
end

