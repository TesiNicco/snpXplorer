from flask import Flask, request, render_template
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io
import math
import pandas as pd
import numpy as np
import base64
import os

app = Flask(__name__)

@app.route("/", methods=["GET", "POST"])

def index():
    if request.method == "POST":
        # read inputs
        # chromosome
        chrom = int(request.form["chrom"])
        # position
        position = int(request.form["position"])
        # gwas
        gwas = request.form.getlist("gwas")
        # window around position
        window = int(request.form["window"])

        # gather data for the plot
        start_pos = position - window
        end_pos = position + window
        df = get_data_plot(gwas = gwas, chrom = chrom, start_pos = start_pos, end_pos = end_pos)

        # plot
        img = scatterplot(df = df, chrom = chrom, start_pos = start_pos, end_pos = end_pos, gwas = gwas)
        # set the plot url for showing on the application
        plot_url = base64.b64encode(img.getvalue()).decode()
        
        # return the html template and the url to the plot
        return render_template("index.html", plot_url=plot_url)
    return render_template("index.html")

# function to get data to plot given the gwas name, chromosome, start and end position
def get_data_plot(gwas, chrom, start_pos, end_pos):
    # set path to data folder
    data_path = '/Users/nicco/Documents/GitHub/snpXplorer/bin/snpxplorer_v4/data'
    # iterate over the gwas to plot
    data_list = []
    for gw in gwas:
        # define command for tabix
        cmd = 'tabix -b 3 -e 3 -s 1 ' + data_path + '/' + gw + '/chr' + str(chrom) + '_' + gw + '.txt.gz ' + str(chrom) + ':' + str(start_pos) + '-' + str(end_pos)
        # extract summary statistics and add gwas name
        sumstats = [[gw] + x.rstrip().split('\t') for x in os.popen(cmd)]
        # append to data_list
        data_list = sumstats if len(data_list) == 0 else data_list + sumstats
    # convert data_list to dataframe
    df = pd.DataFrame(data_list, columns=['Gwas', 'NA', 'Chrom', 'Pvalue', 'Position', 'Rsid', 'Maf', 'Ref', 'Alt'])
    # make position and pvalue numeric
    df["Position"] = pd.to_numeric(df['Position'], errors='coerce')
    df["Pvalue"] = pd.to_numeric(df['Pvalue'], errors='coerce')
    # then take the -log10 of the pvalue
    df['Pvalue'] = -np.log10(df['Pvalue'])
    # then take the position value as megabases
    df["Position"] = df["Position"] / 1000000
    return df

# function to draw the scatterplot given a dataframe and positions to plot
def scatterplot(df, chrom, start_pos, end_pos, gwas):
    # define structure for the plot
    fig, axs = plt.subplots(2, 1) 
    
    # Scatterplot here
    # define colors and map them to datasets
    all_colors = ['red', 'blue', 'orange', 'yellow', 'grey', 'green']
    color_map = {val: all_colors[index] for index, val in enumerate(gwas)}
    colors = df['Gwas'].map(color_map)
    # plot data points
    axs[0].scatter(df['Position'], df['Pvalue'], c=colors)
    # add title
    title = '%s:%s-%s' %(str(chrom), str(start_pos), str(end_pos))
    axs[0].set_title(title)
    # axes names
    axs[0].set_xlabel("")
    axs[0].set_ylabel("-log10(P-value)")
    # disable x-axis values
    axs[0].set_xticks([])
    # Add the legend
    handles = [plt.plot([],[],color=color_map[color], marker="o", ls="", label=color)[0] for color in color_map.keys()]
    axs[0].legend(handles=handles)

    # Gene track here
    axs[1].scatter(df['Position'], df['Pvalue'], c=colors)
    # axes names
    axs[1].set_xticks([])
    axs[1].set_ylabel("Gene track")
    
    # save image
    img = io.BytesIO()
    fig.savefig(img, format="png")
    img.seek(0)
    return img

if __name__ == "__main__":
    app.run()