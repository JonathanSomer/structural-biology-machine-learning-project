import os
import shutil

from Constants import BASE_DATA_PATH

overwrite_existing = False

for root, dirs, files in os.walk(BASE_DATA_PATH):
    for f in files:
        filename, file_extension = os.path.splitext(f)
        fullpath = os.path.join(root, f)
        direct_dir = os.path.basename(os.path.dirname(fullpath))
        if direct_dir  in ('benchmark', 'patch_dock') and 'cache' in filename or file_extension == '.json':
            new_dirpath = os.path.join(os.path.dirname(root), 'cache')
            if not os.path.exists(new_dirpath):
                print("Creating dir {}".format(new_dirpath))
                os.makedirs(new_dirpath)
            new_path = os.path.join(new_dirpath, f)
            if not overwrite_existing and os.path.exists(new_path):
                continue
            print "move {} to {}".format(fullpath, new_path)
            os.rename(fullpath, os.path.join(new_dirpath, f))


