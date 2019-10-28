#script to download file given a url, understand column names, split file in chromosomes and update main server and ui scripts

import sys
import os
import gzip

#function to download file
def download_data(url):
  print("Downloading GWAS dataset...")
  cmd_1 = "mkdir data/tmp0123"
  cmd = "wget %s data/tmp0123/"
  os.system(cmd_1)
  os.system(cmd)
  return ("GWAS dataset downloaded correctly.")
  
#function to understand header
def split_one(MODE, GWAS, TIT):
  #check mode in order to take either the downloaded file or user-downloaded GWAS
  if MODE == "1":
    cmd = "ls data/tmp0123/"
    ls = os.popen(cmd).read()
    ls = ls.rstrip().split()
    f = "data/tmp0123/" + ls[0]
  else:
    f = GWAS
    
  #prepare folder
  cmd = "mkdir data/%s" %(TIT)
  os.system(cmd)

  #create list for pvalue, chromosome and position names to be recognized
  chrom_field = ["chr", "Chr", "CHR", "chrom", "Chrom", "CHROM", "CHROMOSOME"]
  pos_field = ["bp", "BP", "Bp", "pos", "POS", "Pos", "POSITION"]
  p_field = ["p", "P", "pval", "PVAL", "pvalue", "PVALUE", "P_IGAP1", "P_VALUE", "p_value"]

  # parameters: c is general index of lines, h is for header, indexes is to store important column in file, prev_chr is to
  # keep track of chromosomes executed, snp is to keep track of snps per chromosome
  c = 0
  snp = 0
  h = ""
  prev_chr = 0
  indexes = {}
    
  with gzip.open(f) as inpf:
    for line in inpf:
      line = line.rstrip().split()
      
      if (c == 0):
        h = line
        #find indexes for chromosome, position, pvalue
        for i in range(len(h)):
          if h[i] in chrom_field:
            indexes["chr"] = i+1
          elif h[i] in pos_field:
            indexes["pos"] = i+1
          elif h[i] p_field:
            indexes["p"] = i+1
        c = c + 1
      
      else:
        ind_chr, ind_pos, ind_p = indexes["chr"]-1, indexes["pos"]-1, indexes["p"]-1
        chr_n, pos, p = line[ind_chr], line[ind_pos], line[ind_p]
        snp = snp + 1
        
        #check for assigning first chromosome
        if (c == 1):
          prev_chr = chr_n
          c = c + 1
        
        if (prev_chr != chr_n):
          print("### Done with chromsome %s: there were %s SNPs" %(prev_chr, snp))
          prev_chr = chr_n
          snp = 1
          out.close()
        
        #prepare output file
        fname = "data/%s/chr%s_%s.txt" %(TIT, chr_n, TIT)
        out = open(fname, "a")
        
        out.write(str(chr_n))
        out.write(" ")
        out.write(str(pos))
        out.write(" ")
        out.write(str(p))
        out.write("\n")

  return "Done!"

#function to add as repository in server file
def addRepoServer(MODE, TIT):
  #open server script
  srv = open("bin/server.R").readlines()
  
  #new file
  srv2 = open("bin/server.R", "w")
  
  for line in srv:
    if line.startswith("  } #else here to add repositories"):
      c1 = "\t} else if (inp == '%s') {" %(TIT)
      srv2.write(c1)
      srv2.write("\n")
      srv2.write("\t\t")
      c2 = "dat <- fread('../data/%s/chr19_%s.txt', h=T, stringsAsFactors=F)" %(TIT, TIT)
      srv2.write(c2)
      srv2.write("\n")
      srv2.write("\t\t")
      srv2.write("colnames(dat) <- c('chr', 'pos', 'p')")
      srv2.write("\n")
      srv2.write("\t\t")
      srv2.write("dat <- dat[!which(is.na(dat$p)),]")
      srv2.write("\n")
      srv2.write("\t\t")
      srv2.write("chrom <- dat$chr[1]")
      srv2.write("\n")
      srv2.write("\t\t")
      srv2.write("dat$p <- as.numeric(dat$p)")
      srv2.write("\n")
      srv2.write("\t\t")
      srv2.write("dat$'-log10(P-value)' <- -log10(dat$p)")
      srv2.write("\n")
      srv2.write("\t\t")
      c3 = "gwas <- '%s'" %(TIT)
      srv2.write(c3)
      srv2.write("\n")
      c4 = "  } #else here to add repositories"
      srv2.write(c4)
      srv2.write("\n")
      
    else:
      srv2.write(line)
    
  srv2.close()
  return ("Done!")
  
#function to add as repository in ui file
def addRepoUI(MODE, TIT):
  #open server script
  srv = open("bin/ui.R").readlines()
  
  #new file
  srv2 = open("bin/ui.R", "w")
  
  for line in srv:
    if line.startswith("            checkboxGroupInput"):
      line_spl = line.split(",")
      last = line_spl[-3]
      replacement = '"%s" = "%s", "Your input"' %(TIT, TIT)
      repl = last.replace('"Your input"', replacement)
      newL = line_spl[:-3]
      newL.append(repl)
      newL.append(line_spl[-2])
      newL.append(line_spl[-1])
      new_line = ",".join(newL)
      
      srv2.write(new_line)
      srv2.write("\n")
    else:
      srv2.write(line)
    
  srv2.close()
  return ("Done!")

#function to clean data if already present (take only 3 columns to save some space and speed up)
def slimGWAS(GWAS, MODE, TIT):
  #prepare folder
  cmd = "mkdir data/%s" %(TIT)
  os.system(cmd)
  
  #create list for pvalue, chromosome and position names to be recognized
  chrom_field = ["chr", "Chr", "CHR", "chrom", "Chrom", "CHROM"]
  pos_field = ["bp", "BP", "Bp", "pos", "POS", "Pos"]
  p_field = ["p", "P", "pval", "PVAL", "pvalue", "PVALUE", "P_IGAP1", "P_VALUE", "p_value"]

  #main loop on chromosomes
  for chrom in range(1, 23):
    #uder update
    print("\t## Working on chromosome %s"%(chrom))
    
    #this assumes the data are already divided in chromosomes, but there are redundant columns that slow down the program
    #then look on the chromosomes -- GWAS is the chromosome1 file for the GWAS -- also assume the pattern chrXX_blabla (split by chr and _)
    fname = GWAS.split("chr")
    fname2 = fname[-1]
    fname2 = fname2.split("_")
    fname3 = "_".join(fname2[1:])
    path = fname[0] + "chr" + str(chrom) + "_" + fname3

    #check if zipped
    flag_zip = False
    if fname3.endswith("zip") or fname3.endswith("gz"):
      flag_zip = True
    
    # parameters: c is general index of lines, h is for header, indexes is to store important column in file, prev_chr is to
    # keep track of chromosomes executed, snp is to keep track of snps per chromosome
    c = 0
    snp = 0
    h = ""
    prev_chr = 0
    indexes = {}
    
    print(path)
    
    #open file (take zip into account)
    if flag_zip == True:
      with gzip.open(path) as inpf:
        for line in inpf:
          line = line.rstrip().split()

          if (c == 0):
            h = line
            print(h)
            #find indexes for chromosome, position, pvalue
            for i in range(len(h)):
              if h[i] in chrom_field:
                indexes["chr"] = i+1
              elif h[i] in pos_field:
                indexes["pos"] = i+1
              elif h[i] in p_field:
                indexes["p"] = i+1
            c = c + 1
      
          else:
            ind_chr, ind_pos, ind_p = indexes["chr"]-1, indexes["pos"]-1, indexes["p"]-1
            chr_n, pos, p = line[ind_chr], line[ind_pos], line[ind_p]
            snp = snp + 1
        
            #check for assigning first chromosome
            if (c == 1):
              prev_chr = chr_n
              c = c + 1
        
            if (prev_chr != chr_n):
              print("### Done with chromsome %s: there were %s SNPs" %(prev_chr, snp))
              prev_chr = chr_n
              snp = 1
              out.close()
          
            #prepare output file
            fname = "data/%s/chr%s_%s.txt" %(TIT, chr_n, TIT)
            out = open(fname, "a")
        
            out.write(str(chr_n))
            out.write(" ")
            out.write(str(pos))
            out.write(" ")
            out.write(str(p))
            out.write("\n")
          
    else:
      with open(path) as inpf:
        for line in inpf:
          line = line.rstrip().split()
          
          if (c == 0):
            h = line
            #find indexes for chromosome, position, pvalue
            for i in range(len(h)):
              if h[i] in chrom_field:
                indexes["chr"] = i+1
              elif h[i] in pos_field:
                indexes["pos"] = i+1
              elif h[i] in p_field:
                indexes["p"] = i+1
            c = c + 1
      
          else:
            ind_chr, ind_pos, ind_p = indexes["chr"]-1, indexes["pos"]-1, indexes["p"]-1
            chr_n, pos, p = line[ind_chr], line[ind_pos], line[ind_p]
            snp = snp + 1
        
            #check for assigning first chromosome
            if (c == 1):
              prev_chr = chr_n
              c = c + 1
        
            if (prev_chr != chr_n):
              print("\t## Done with chromsome %s: there were %s SNPs" %(prev_chr, snp))
              prev_chr = chr_n
              snp = 1
              out.close()
          
            #prepare output file
            fname = "data/%s/chr%s_%s.txt" %(TIT, chr_n, TIT)
            out = open(fname, "a")
        
            out.write(str(chr_n))
            out.write(" ")
            out.write(str(pos))
            out.write(" ")
            out.write(str(p))
            out.write("\n")
  return "Done!"          
          
#################################
########### MAIN ################
#################################

#running mode: 1=download+split+parse; 2=split+parse; 3=parse
MODE = sys.argv[1]

#healp mode
if MODE == "--help":
  print("##### ADD NEW GWAS #####")
  print("### ARGUMENTS:")
  print("### arg1: MODE (--download-split-parse=download+split+add; --split-add=split+add; --add=add; --slim-add=slim+add; --help=this_message)")
  print("### Supplementary arguments:")
  print("### \tMODE --download-split-add: arg2 --> URL (The url of the GWAS to download)")
  print("### \tMODE --download-split-add: arg3 --> TIT (This is the title for the GWAS and folders)")
  print("#############################################################################")
  print("### \tMODE --split-add: arg2 --> GWAS (This is the name of the GWAS dataset)")
  print("### \tMODE --split-add: arg3 --> TIT (This is the name of the folder and GWAS)")
  print("#############################################################################")
  print("### \tMODE --add: arg2 --> TIT (This is the name of the folder and GWAS)")
  print("#############################################################################")
  print("### \tMODE --slim-add: arg2 --> GWAS (This is the name of the chromosome 1 of the GWAS dataset)")
  print("### \tMODE --slim-add: arg3 --> TIT (This is the name of the folder and GWAS)")
  print("#############################################################################\n\n")
  
#running mode download+split+add
elif MODE == "--download-split-parse":
  #get argument 2
  URL = sys.argv[2]
  
  #get title as argument 3
  TIT = sys.argv[3]

  #create fake GWAS variable
  GWAS = ""
  
  #update
  print("## ADD NEW GWAS ##")
  print("##################")
  print("##################")
  print("## You choose:")
  print("## Mode = --download-split-add")
  print("## URL = %s" %(URL))
  print("## GWAS name = %s (name of folder in data/ and each file)" %(TIT))
  print("##################")
  print("##################\n\n")

  #download gwas dataset
  print (download_data(url))

  #understand structure of gwas file
  print split_one(MODE, GWAS, TIT)
  
  #clean tmp folder
  cmd = "rm -rf data/tmp0123"
  os.system(cmd)

  #need to add as repository -- server file
  print addRepoServer(MODE, TIT)
  
  #need to add as -- ui file
  print addRepoUI(MODE, TIT)
  
#running mode split+add
elif MODE == "--split-add":
  
  #argument 1 is GWAS path
  GWAS = sys.argv[2]
  
  #get title as argument 3
  TIT = sys.argv[3]

  #update
  print("## ADD NEW GWAS ##")
  print("##################")
  print("##################")
  print("## You choose:")
  print("## Mode = --split-add")
  print("## GWAS file = %s (name of GWAS dataset)" %(GWAS))
  print("## GWAS name = %s (name of folder in data/ and each file)" %(TIT))
  print("##################")
  print("##################\n\n")
  
  #split in file per chromosome after identifying the header
  print split_one(MODE, GWAS, TIT)

  #need to add as repository -- server file
  print addRepoServer(MODE, TIT)
  
  #need to add as -- ui file
  print addRepoUI(MODE, TIT)

#running mode slim+add
elif MODE == "--slim-add":
  #path to chromosome 1 of GWAS is the first argument
  GWAS = sys.argv[2]

  #title is the second argument
  TIT = sys.argv[3]
  
  #update
  print("## ADD NEW GWAS ##")
  print("##################")
  print("##################")
  print("## You choose:")
  print("## Mode = --slim ")
  print("## GWAS file = %s (name of GWAS dataset)" %(GWAS))
  print("## GWAS name = %s (name of folder in data/ and each file)" %(TIT))
  print("##################")
  print("##################\n\n")

  #run function to slim data
  print(slimGWAS(GWAS, MODE, TIT))

  #need to add as repository -- server file
  print addRepoServer(MODE, TIT)
  
  #need to add as -- ui file
  print addRepoUI(MODE, TIT)

#running mode add only 
elif MODE == "--add":
  #title is the second argument
  TIT = sys.argv[2]

  #update
  print("## ADD NEW GWAS ##")
  print("##################")
  print("##################")
  print("## You choose:")
  print("## Mode = --add ")
  print("## GWAS name = %s (name of folder in data/ and each file)" %(TIT))
  print("##################")
  print("##################\n\n")
  
  #need to add as repository -- server file
  print addRepoServer(MODE, TIT)
  
  #need to add as -- ui file
  print addRepoUI(MODE, TIT)
  
  
  
  
  
  
  
  
  
  
