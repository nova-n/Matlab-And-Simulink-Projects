#include <iostream>
#include <vector> //so can add elements to arrays at runtime
#include <cmath>

using namespace std;
int desiredDigits = 3;

vector<int> primeFinder(int digitCount){
    int maxNum = pow(10,digitCount) - 1;
    //for example, if want 3 digits, then get 10^3 -1 = 1000 -1 = 999, the largest 3 digit numbner
    vector<int> primes = {2,3};

    //cout << "maxNum: "<<maxNum << " \n";
    int ii; //need to declare this outside the loop so that stuff outside the loop can access it

    //sizeof(primes)/sizeof(primes[0]) is the length of the array, primes with plain c++ arrays
    //primes.size() is how it is done with vectors in c++
    //will start collecting at the end of the primes list + 1, so in this case, at 3+1=4
    for(int i = primes[primes.size() - 1] +1; i<=maxNum ; i++){
        //cout << "i = " << i << " \n";
        for(ii = 0;ii<primes.size();ii++){
            if(i % primes[ii] == 0){
                break;
                //modulo to check for remainders of dividing i by primes(ii)
                //if remainder is 0, then that number is divisible, so break
                //loop, and skip this number
            }
        }
        //cout << "checked, ii = " << ii << " , length of primes: " << primes.size() << " \n";
        //if it reached the end of the list without breaking loop (for ii), then
        //that number was not found to be divisible by anything, so its prime
        //this works because you keep adding more and more numbers to the list,
        //so bigger numbers have more factors to check against
        if(ii == primes.size()){
            primes.push_back(i);
        }
    }
    return primes;
}

int main(){
    int sumOfPrimes=0;
    vector<int> foundPrimes = primeFinder(desiredDigits);
    cout << "Primes Up To " << desiredDigits << " Digits: \n";
    for (const int &i : foundPrimes) {
        cout << i << " , ";
        sumOfPrimes += i;
    }
    cout << "\n";
    cout << "Sum Of Primes Up To " << desiredDigits << " Digits: "<< sumOfPrimes <<"\n";
    return 0;
}

