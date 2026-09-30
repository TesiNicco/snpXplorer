# Script to compress and check PRS (Polygenic Risk Scores) data

# conda activate nlpy

# Libraries
import os
import pandas as pd

# Data path
data_path_euro = 'OpenGWAS_sumstats/European/GWASum'
data_path_african = 'OpenGWAS_sumstats/African/GWASum'

# List data in both directories
list_cmd_euro = f"ls {data_path_euro}"
list_cmd_african = f"ls {data_path_african}"

# Execute the list commands
list_euro = [x for x in os.popen(list_cmd_euro).read().splitlines() if x != "reformatted"]
len(list_euro)
list_african = [x for x in os.popen(list_cmd_african).read().splitlines() if x != "reformatted"]
len(list_african)

# Then take the list of files after GWASum reformatted
list_euro_reformatted = [x.replace('_clumped_reformatted.txt', '') for x in os.popen(f"ls {data_path_euro}/reformatted").read().splitlines()]
list_african_reformatted = [x.replace('_clumped_reformatted.txt', '') for x in os.popen(f"ls {data_path_african}/reformatted").read().splitlines()]
len(list_euro_reformatted)
len(list_african_reformatted)

# See what files did not get reformatted
not_reformatted_euro = [x for x in list_euro if x not in list_euro_reformatted]
not_reformatted_african = [x for x in list_african if x not in list_african_reformatted]
len(not_reformatted_euro)
len(not_reformatted_african)    # 2 -- peace, will skip

# Iterate over the files that did get reformatted and compress them
counter = 0
for file in list_euro_reformatted:
    counter += 1
    print(f"[{counter}/{len(list_euro_reformatted)}] Compressing {file}")
    full_path = os.path.join(data_path_euro, "reformatted", f"{file}_clumped_reformatted.txt")
    if os.path.exists(full_path):
        os.system(f"gzip {full_path}")
counter = 0
for file in list_african_reformatted:
    counter += 1
    print(f"[{counter}/{len(list_african_reformatted)}] Compressing {file}")
    full_path = os.path.join(data_path_african, "reformatted", f"{file}_clumped_reformatted.txt")
    if os.path.exists(full_path):
        os.system(f"gzip {full_path}")

# Next step is to combine all traits and ids from European and African together, and add the ancestry
# Define data path
african_traits_path = "/project/holstegelab/Share/nicco/haplotypeBlocks/R2cluster/OpenGWAS_sumstats/African/trait_names_and_ids.txt"
european_traits_path = "/project/holstegelab/Share/nicco/haplotypeBlocks/R2cluster/OpenGWAS_sumstats/European/trait_names_and_ids.txt"

# Read the trait names and IDs from both European and African files and add the ancestry information
african_traits_df = pd.read_csv(african_traits_path, sep='\t')
african_traits_df['ancestry'] = 'African'
european_traits_df = pd.read_csv(european_traits_path, sep='\t')
european_traits_df['ancestry'] = 'European'
# Combine
combined_traits_df = pd.concat([african_traits_df, european_traits_df], ignore_index=True)
combined_traits_df.to_csv("/project/holstegelab/Share/nicco/haplotypeBlocks/R2cluster/OpenGWAS_sumstats/combined_traits_and_ids.txt", sep='\t', index=False)