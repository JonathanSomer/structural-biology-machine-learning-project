from Constants import *
import re
import json


def main():
	for complex_id in get_all_training_complexes():
		rank = 1
		map_rank_to_score_components = {}
		
		print("STARTED BUILDING {}".format(complex_id))
		with open(get_patchdock_complex_score_file_path(complex_id), "r") as f:
		    for line in f:
		        pattern = re.compile("^\s*" + str(rank) + "\s\|.+")
		        if pattern.match(line):
		            components = [x.strip() for x in line.split('|')]
		            # rank | score | pen.  | Area    | as1   | as2   | as12  | ACE     | hydroph | Energy  |cluster| dist. || Ligand Transformation
		            indexes = [1, 2, 3, 7]
		            components = [float(components[i]) for i in indexes]
		            map_rank_to_score_components[rank] = components
		            rank += 1
		            if rank > 1000:
		            	break

		# print(map_rank_to_score_components)
		with open(get_patchdock_complex_score_json_file_path(complex_id), "w") as f:
			json.dump(map_rank_to_score_components, f)

