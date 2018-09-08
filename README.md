# Structural Biology Machine Learning Project

## Links:
1. [Tasks Doc](https://docs.google.com/document/d/1qvxntDFh5F84yEO1LBVDT8cNZrA7_KstXjd-2RdRQLM/edit?ts=5b67f9e8)
2. [Zdock Benchmark](https://zlab.umassmed.edu/benchmark)
3. [RaptorX](http://raptorx.uchicago.edu/ComplexContact/)
4. [data folder - must download and place in home in order to run stuff](https://drive.google.com/open?id=1XHeE9714OMDclEfTrm1s_RUv1-ve3Jo1)


## Guide for the 'working with nova' scripts:
1. First step is necassary only for working on a new server, run:

	`scp -r path/to/PatchDock/directory/in/repo <USERNAME>@nova.cs.tau.ac.il:`

2. Transfer all benchmark complexes to server and run patch dock on them:

	```python
	cd nova/
	python3 run_patch_dock_on_entire_benchmark.py -u <USERNAME> -p <PASSWORD>
	```

3. Now in order to work with transformations locally, insert the caching logic into the relevant location in transformation_getter.py and run:

	```python
	cd nova/
 	python3 transformation_getter.py -c 2A9K -n 1000 -u <USERNAME> -p <PASSWORD>
	```

	note: you can use the -k flag in order to keep the pdb files when done
