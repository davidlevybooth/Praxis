#!/bin/bash

read -p "Enter a name for your pipeline: " name

conda env create --name $name --file envs/environment.yaml

if [ $(grep 'Praxis' ~/.bashrc | wc -l) -eq 0 ]
  then
    cat Praxis_function.txt >> ~/.bashrc
    echo "You must start a new terminal session for the Praxis command to work."
fi

