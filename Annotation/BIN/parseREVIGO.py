#!/usr/bin/python
# PACKAGES
import re
import werkzeug
werkzeug.cached_property = werkzeug.utils.cached_property
from robobrowser import RoboBrowser
import sys

# MAIN
## read random number to take revigon input
random_num = sys.argv[3]
#random_num = 4167

## then define input file
inpf = "RESULTS_" + str(random_num) + "/revigo_inp.txt"

## read inputs go and pvalues
d = open(inpf).read()

## set connections
br = RoboBrowser(parser="lxml")
#br = RoboBrowser(parser="html5lib")
br.open("http://revigo.irb.hr/")

## manage parameters
clus_sim = float(sys.argv[1])
if clus_sim == 0.4:
    clus_sim = '0.40'
elif clus_sim == 0.5:
    clus_sim = '0.50'
elif clus_sim == 0.7:
    clus_sim = '0.70'
#clus_sum = '0.4'

sim_meas = sys.argv[2]
if sim_meas == "Lin":
    sim_meas = 'Lin'
if sim_meas == "SimRel":
    sim_meas = 'SIMREL'
#clus_sim = 'Lin'

## fill in my data
form = br.get_form()
#form["ctl00$MasterContent$txtGOInput"] = d

form["goList"].value = d
form["cutoff"].value = clus_sim
form["goSizes"].value = '9606'
form["measure"].value = sim_meas

## submit request
br.submit_form(form)

## get link to download stuff
download_csv_link = br.find("a", href=re.compile("export.jsp"))
br.follow_link(download_csv_link)
csv_content = br.response.content.decode("utf-8")

## define output file name
fout_name = "RESULTS_" + str(random_num) + "/revigo_out.csv"

## write output
fout = open(fout_name, "w")
fout.write(csv_content)
fout.close()



