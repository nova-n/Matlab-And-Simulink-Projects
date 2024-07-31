#!/bin/bash

echo "ALLAHU AKBAR BITCH"

#array of naming extensions
Re=('400' '1000' '3200') #no spaces!
meshes=('20x20' '40x40' '60x60')


for i in {0..2} # ${#Re[@]}
do
    for ii in {0..2} # ${#meshes[@]}
    do 
        fileName="cavity_Re${Re[i]}_mesh${meshes[ii]}"
        #x = "cavityFiles/physicalProperties" #Doesn't actually modify the number itself, would have to do that on your own
        echo $fileName
        cp -a cavity/. $fileName
        #cp $
        #cp

        # ----- run simulation
        cd $fileName
        ./run.sh #need to copy this run file over from outside
    done
done


