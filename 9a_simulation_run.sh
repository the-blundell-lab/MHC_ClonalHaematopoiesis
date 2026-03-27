#!/bin/bash

# this script allows to run the python code in a for loop (for the same mu and N)
# with different values of the s (fitness) coefficient in 0.5% intervals 
# results saved to folder msc_thesis/results/simulation/run_across_s

# RUN: 
# cd ~/Desktop/msc_thesis/task1_predict_binding_to_HLA/scripts02
# bash 9a_simulation_run.sh

# Specify path of the folder to save results to 
save_folder="/Users/barbarawalkowiak/Desktop/msc_thesis/results/simulation"

# Define start (min s = 10%), stop (max s = 20%) and difference in s (0.25%)
START_s=0.1
END_s=0.2
STEP_s=0.0025

# for mu in 1e-9 2e-9 5e-9 1e-8 2e-8 5e-8
for mu in 5e-7

do

    current=$START_s
    
    echo "mu: $mu"

    # Loop through the numerical range
    while (( $(echo "$current <= $END_s" | bc -l) )); do

        echo "s: $current"

        # avoid running this if the simulation already exists 
        folder_name=simulation_"${current}"_"${mu}"

        echo "$folder_name"

        if [ ! -d "$save_folder/$folder_name" ]; then

            echo "Running simulation with fitness coefficient: $current and mutation rate: $mu"
            python3 9a_simulation_code.py --name simulation_"$current"_"$mu" --mu "$mu" --s "$current" --N 100000 --depth 70 --lifespan 74

            # get the new current value (current + step)
            current=$(echo "$current + $STEP_s" | bc -l)

        else

            echo "combination found. not running this simulation again"
            current=$(echo "$current + $STEP_s" | bc -l)
            echo "$current"

        fi 

    done

done