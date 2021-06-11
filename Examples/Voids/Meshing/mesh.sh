#!/bin/bash

# Clear working directory
rm -r mesh_scratch
mkdir mesh_scratch
cd mesh_scratch

# Copy compiled files
cp ../../../../crackMesher .
cp ../../../../stl_binary2ascii.py .
cp ../../../shrinkVol.py .

cp ../config.xml .
cp ../../VoidList .
cp ../../holes .

# Copy input stl files
cp -R ../../STL ./stls

# Run the code
./crackMesher config.xml

# Scale the model to the correct size
python3 shrinkVol.py
cd ..