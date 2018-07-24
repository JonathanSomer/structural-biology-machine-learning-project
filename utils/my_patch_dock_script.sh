cd /specific/a/home/cc/students/csguests/somer/PatchDock
./buildParams.pl /specific/a/home/cc/students/csguests/somer/patch_dock_results/benchmark/structures/1AHW_r_u.pdb /specific/a/home/cc/students/csguests/somer/patch_dock_results/benchmark/structures/1AHW_l_u.pdb

cd /specific/a/home/cc/students/csguests/somer/patch_dock_results/data
mkdir 1AHW
cd 1AHW
mkdir patch_dock_ranking
mkdir benchmark

cd /specific/a/home/cc/students/csguests/somer/PatchDock
patch_dock.Linux params.txt /specific/a/home/cc/students/csguests/somer/patch_dock_results/data/1AHW/patch_dock_ranking/1AHW

transOutput.pl /specific/a/home/cc/students/csguests/somer/patch_dock_results/data/1AHW/patch_dock_ranking/1AHW 1 1000
cd /specific/a/home/cc/students/csguests/somer/patch_dock_results/data/1AHW/patch_dock_ranking/

scp -r somer@nova.cs.tau.ac.il:/specific/a/home/cc/students/csguests/somer/patch_dock_results/data/1AHW /Users/jonathan.somer/School/structural-biology-machine-learning-project/data/
