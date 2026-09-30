# Script to gather SNP associations from OpenGWAS based on population and sample size criteria

# Libraries
import pandas as pd
from ieugwaspy import query

# Define population of interest and the relative output path
population_of_interest = "African"
if population_of_interest == "European":
    output_folder = "OpenGWAS_sumstats"
elif population_of_interest == "African":
    population_of_interest = ["African American or Afro-Caribbean", "Sub-Saharan African", "African unspecified"]
    output_folder = "OpenGWAS_sumstats/African"

# Define path to folder where to save sumstats of interest

# List all studies, then (optionally) filter by ancestry and sample size
studies = pd.DataFrame(query.gwasinfo())  # returns ALL studies' metadata

# Transpose so each row = one study
studies = studies.T  # <- transpose

# Reset index and rename
studies.reset_index(inplace=True)

# Check population and sample size
studies["population"].value_counts()
studies["sample_size"] = pd.to_numeric(studies["sample_size"], errors="coerce")

# Restrict to studies based on ancestry -- for European, add also sample size > 10000 - for African, sample size > 5000
studies_filtered = studies[studies["population"].isin(population_of_interest)]
if population_of_interest == "European":
    studies_filtered = studies_filtered[studies_filtered["sample_size"] > 10000]
elif population_of_interest == ["African American or Afro-Caribbean", "Sub-Saharan African", "African unspecified"]:
    studies_filtered = studies_filtered[studies_filtered["sample_size"] > 5000]
studies_filtered["sample_size"].describe()

# Drop duplicates based on id
studies_filtered = studies_filtered.drop_duplicates(subset=["id"])
list_ids = studies_filtered["id"].tolist()

# Extract SNPs with association p-value < 5e-5 across all filtered studies
counter = 0
for sid in list_ids:
    counter += 1
    try:
        hits = query.tophits(id=[sid], pval=5e-5, clump=0)  # clump=0 = all hits, not LD-pruned
        if hits:
            df_hits = pd.DataFrame(hits)
            df_hits["study_id"] = sid
            # write to file
            df_hits.to_csv(f"{output_folder}/{sid}_hits_p5e-5.txt", sep="\t", index=False)
            print(f"[{counter}/{len(list_ids)}] {sid}: {len(df_hits)} hits")
    except Exception as e:
        print(f"Skipping {sid}: {e}")
