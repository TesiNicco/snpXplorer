# script to liftover new summary statistics assuming hg38 data per chromosome is there with the right columns

# libraries
from liftover import get_lifter
import os
import pandas as pd

# read files to be lifted
directory_files_to_lift = '/root/snpXplorer/data/Bellenguez'
flist = [x.rstrip() for x in os.popen('ls %s/*gz' %(directory_files_to_lift)).readlines() if 'chr' in x]

# iterate across chromosomes
for f in flist:
    print('** working on %s' %(f))
    # read data into dataframe
    df = pd.read_csv(f, sep='\t')
    # liftover
    converter = get_lifter('hg38', 'hg19')
    # liftover
    lifted_positions = []
    for index, row in df.iterrows():
        pos_interest = converter[df['CHR'][index]][df['POS'][index]]
        if len(pos_interest) >0:
            lifted_positions.append(pos_interest[0][1])
        else:
            lifted_positions.append('')
    # then add the list as additional column
    df['POS_HG37'] = lifted_positions
    # create a new df
    columns_to_select = ['POS_HG37', 'CHR', 'P', 'POS', 'RSID', 'MAF', 'REF', 'ALT']
    new_df = df[columns_to_select]
    new_df.columns = ['POS', 'CHR', 'P', 'POS_HG38', 'RSID', 'MAF', 'REF', 'ALT']
    # exclude empty --> liftover failed and make numeric
    new_df = new_df[new_df['POS'] != '']
    new_df['POS'] = pd.to_numeric(new_df['POS'], errors='coerce')
    # sort by position
    sorted_df = new_df.sort_values(by='POS', ascending=True)
    # finally write data
    outname = f.replace('_hg38', '').replace('.gz', '')
    sorted_df.to_csv(outname, sep='\t', index=False)
    # compress with bgzip
    os.system('bgzip %s' %(outname))
    # tabix
    os.system('tabix -S 1 -s 2 -b 1 -e -1 %s.gz' %(outname))

        

