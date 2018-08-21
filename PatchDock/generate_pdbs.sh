#!/bin/bash

echo "generating pdbs..."
# complexes=( 1MLC 1VFB 1WEJ 2FD6 2I25 2VIS 2VXT 2W9E 3EOA 3HMX 3MXW 3RVW 4DN4 4FQI 4G6J 4G6M 4GXU 1BJ1 1FSK 1I9R 1IQD 3L5W 3V6Z 1CGI 1IJK 1JIW 1KKL 1M10 1NW9 1R6Q 1ZM4 2NZ8 2Z0E 1ACB 1F6M 1FQ1 1JK9 )
complexes=( 1AVX 1K4C 1KXQ 1NCA 1NSN 1QFW 2JEL )
for complex in "${complexes[@]}"
do
	echo "generating complex for: $complex"
	date -u
	./transOutput.pl ../data/${complex}/patch_dock/${complex}.patch_dock_output 1 200
done
