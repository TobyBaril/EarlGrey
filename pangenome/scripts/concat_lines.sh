#!/bin/bash

# Read the first three lines from each input file and concatenate them
cat "$1" | tr '\n' ' ' > "$2"

#Create a random file
echo " This is a random file." > "${2}_random.txt"