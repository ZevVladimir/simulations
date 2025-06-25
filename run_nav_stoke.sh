#!/bin/bash
rm -r ~/simulations/nav_stoke_output/*

make clean
make nav_stoke

mkdir ~/.screen && chmod 700 ~/.screen
export SCREENDIR=$HOME/.screen

screen -dmS Navier_Stoke_Run 
./nav_stoke

/home/zvladimi/.pyenv/versions/mlois_env/bin/python /home/zvladimi/simulations/plot_nav_stoke.py