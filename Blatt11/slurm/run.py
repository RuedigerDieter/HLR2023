
import os
# assign directory
directory = '.'

# iterate over files in
# that directory
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    # checking if it is a file
    if os.path.isfile(f):
        if filename.endswith(".out"):
            print("./gather_results.sh " + filename)
            os.system("./gather_results.sh " + filename)
            #os.system("sbatch slurm/" + filename)