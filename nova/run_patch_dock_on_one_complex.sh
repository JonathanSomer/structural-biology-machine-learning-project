#!/bin/bash
echo "STARTING PATCH DOCK SCRIPT"
complex=$1
echo "STARTING TO WORK ON: $complex"
date -u
./buildParams.pl ../data/${complex}/benchmark/${complex}_r_u.pdb ../data/${complex}/benchmark/${complex}_l_u.pdb 4.0 EI
./patch_dock.Linux params.txt ${complex}.patch_dock_output
mkdir ../data/${complex}/patch_dock/
mv ${complex}.patch_dock_output ../data/${complex}/patch_dock/
