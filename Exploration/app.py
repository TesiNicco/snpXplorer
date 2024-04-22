from flask import Flask, request, render_template, send_file, session, make_response, redirect, url_for, jsonify
from flask_sqlalchemy import SQLAlchemy
from flask_session import Session
from werkzeug.security import generate_password_hash, check_password_hash
import redis
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io
from io import StringIO
import subprocess
import math
import pandas as pd
import re
import numpy as np
import webcolors
import random
import jinja2
import matplotlib.colors as mcolors
import base64
import random
import statsmodels.api as sm
from liftover import get_lifter
import os
# differences between local and server
# local
#from pandasgwas.get_SNPs import get_variants_by_variant_id
#data_path = '/Users/nicco/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/GitHub/snpXplorer/bin/snpxplorer_v4/data'
# server
data_path = '/Data'
from pandasgwas import get_variants_by_variant_id

# Initialize the App
app = Flask(__name__)
# Add configuration for the SQLAlchemy database (for usernames and passwords)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///db.sqlite'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)
# Add configuration for the Session to save to redis
app.config['SESSION_TYPE'] = 'redis'
app.config['SESSION_REDIS'] = redis.from_url('redis://localhost:6379')
app.config['SESSION_USE_SIGNER'] = True
app.config['SECRET_KEY'] = 'secret'
Session(app)

# Create a user class that contains id, username and password
class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(20), unique=True, nullable=False)
    password = db.Column(db.String(64), nullable=False)
    email = db.Column(db.String(120), unique=True, nullable=True)

# Create database
with app.app_context():
    db.create_all()

# Function to clear the database -- to do so, connext to http://localhost:5000/clear_database
@app.route('/clear_database')
def clear_database():
    db.reflect() # Reflect all tables to Metadata
    db.drop_all() # Drop all those tables
    return "Database Cleared!"

# Registration route
@app.route('/register', methods=['GET', 'POST'])
def register():
    success = None  # set default value for success message
    if request.method == 'POST':
        # getting form data
        username = request.form.get('username')
        password = request.form.get('password')
        email = request.form.get('email')
        confirm_password = request.form.get('confirm_password')

        # validating the form data
        if len(username) < 5:
            error = 'Username must have at least 5 characters.'
            return render_template('user.html', error=error)
        elif len(password) < 8:
            error = 'Password must have at least 8 characters.'
            return render_template('user.html', error=error)
        elif password != confirm_password:
            error = 'Passwords do not match.'
            return render_template('user.html', error=error)

        # hashing and salting the password
        hashed_password = generate_password_hash(password, method='sha256')

        # creating a new user object and adding to the db
        new_user = User(username=username, password=hashed_password, email=email)
        try:
            db.session.add(new_user)
            db.session.commit()
            success = "You have successfully registered. You can Log In now."
        except:
            db.session.rollback()
            error = "Something went wrong during registration. Please try again."
            return render_template('user.html', error_register=error)

        # return success message and render user.html again with success message
        return render_template('user.html', success_register=success)

    return render_template('user.html')

# Login route
@app.route('/login', methods=['GET', 'POST'])
def login():
    # Check if the use already logged in
    if 'user_id' in session:
        return redirect(url_for('dashboard'))
    
    success = None
    error = None
    
    if request.method == 'POST':
        # getting form data
        username = request.form.get('username')
        password = request.form.get('password')

        # querying user from the db
        user = User.query.filter_by(username=username).first()

        # checking if user exists and if password is correct
        if user and check_password_hash(user.password, password):
            success = "Welcome back %s!" %(username)
            # Save the user's ID in the session
            session['user_id'] = user.id
            session['success_login'] = success
            return redirect(url_for('dashboard'))
        else:
            error = "Invalid username or password"

    return render_template('user.html', error_login=error, success_login=success)

# Redirect to user's dashboard
@app.route('/dashboard')
def dashboard():
   return render_template('dashboard.html')

# Index tab
@app.route("/")
def index():
    return render_template("index.html")

# Exploration tab
@app.route('/exploration/', methods=["GET", "POST"])
def exploration():
    # Check if there are new inputs or it is still the previous run
    if request.method == "POST":
        # Read inputs
        # reference genome
        refGen = request.form["refGenome"]
        # browsing option
        browse = request.form["browse"]
        # gwas
        gwas = request.form.getlist('gwas_data')
        # window around position
        window = int(request.form["window"])
        # radio button of sv
        selected_sv_source = request.form.getlist("sv_source")
        # recombination rates
        recomb = request.form["recomb"]
        # show exons
        exons = request.form["exons"]
        # plot type
        plotype = request.form["plotype"]
        # QTL tissues
        qtl_tissues = request.form.getlist('QTLtissuesExplo')

        # gather data for the plot
        # chromosome and positions based on the browsing option
        chrom, start_pos, end_pos, browse_type = readBrowseOption(data_path, browse, window, refGen)
        # get GWAS data to plot
        df = get_data_plot(data_path = data_path, gwas = gwas, chrom = chrom, start_pos = start_pos, end_pos = end_pos, refGen = refGen)
        # get genes to plot
        genes = extract_genes(data_path, chrom, start_pos, end_pos, refGen)
        # get svs to plot
        svs, svs_df = extract_sv(data_path, chrom, start_pos, end_pos, selected_sv_source, refGen)
        # if requested, extract recombination rates
        recomb_data = extract_recomb(data_path, chrom, start_pos, end_pos, refGen) if recomb == 'Yes' else "None"
        # gather data for GTEx
        gtex_df = get_gtex(data_path, genes, refGen)

        # Plot
        img = scatterplot(df = df, chrom = chrom, start_pos = start_pos, end_pos = end_pos, gwas = gwas, genes = genes, svs = svs, browse_type = browse_type, refGen = refGen, recomb_data = recomb_data, exons = exons, browse = browse, plotype = plotype)
        # set the plot url for showing on the application
        plot_url = base64.b64encode(img.getvalue()).decode()
        # then generate the GTEx plot
        img_gtex = gtex_heatmap(gtex_df)
        # set the plot url for showing on the application
        plot_gtex = base64.b64encode(img_gtex.getvalue()).decode()
        
        # Tables
        # prepare the table of SNP associations: sort, subset of columns, add locus and alleles columns, and round Pvalue
        df_sorted = df.sort_values(by='Pvalue', ascending=False)
        df_sorted['Locus'] = df_sorted.apply(lambda x: str(x['Chrom']) + ':' + str(x['Position']), axis=1)
        df_sorted['Alleles'] = df_sorted.apply(lambda x: str(x['Ref']) + '/' + str(x['Alt']), axis=1)
        df_sorted = df_sorted[['Locus', 'Rsid', 'Gwas', 'Alleles', 'Pvalue']]
        df_sorted['Pvalue'] = df_sorted['Pvalue'].round(2)
        snps_table = df_sorted.to_dict(orient='records')
        # prepare the table of the SVs
        svs_df['Region'] = svs_df.apply(lambda x: str(x['chrom']) + ':' + str(x['start']) + '-' + str(x['end']), axis=1)
        svs_df = svs_df[['Region', 'diff', 'type']]
        svs_table = svs_df.to_dict(orient='records')
        # get the table of the gwas catalog
        gwas_table, gwascat_df = extract_gwascatalog(data_path, chrom, start_pos, end_pos, refGen)

        # Store data to enable download of the tables
        session['gtex_df'] = gtex_df.to_csv(index=False)
        session['df'] = df.to_csv(index=False)
        session['svs_df'] = svs_df.to_csv(index=False)
        session['gwascat_df'] = gwascat_df.to_csv(index=False)
        session['chrom'] = chrom
        session['start_pos'] = start_pos
        session['end_pos'] = end_pos
        session['gwas'] = gwas
        session['genes'] = genes
        session['svs'] = svs
        session['browse_type'] = browse_type
        session['refGen'] = refGen
        session['recomb_data'] = recomb_data.to_csv(index=False) if recomb == "Yes" else "None"
        session['exons'] = exons
        session['browse'] = browse
        session['plotype'] = plotype

        # return the html template and the url to the plot
        return render_template("exploration.html", plot_url=plot_url, plot_gtex=plot_gtex, table_snps=snps_table, table_svs=svs_table, table_gwas=gwas_table, browse_value=browse, gwas=gwas)
    elif 'df' in session:
        chrom = session['chrom']
        start_pos = session['start_pos']
        end_pos = session['end_pos']
        gwas = session['gwas']
        genes = session['genes']
        svs = session['svs']
        browse_type = session['browse_type']
        refGen = session['refGen']
        recomb_data = session['recomb_data']
        exons = session['exons']
        browse = session['browse']
        plotype = session['plotype']
        df = pd.read_csv(StringIO(session.get('df')))
        gtex_df = pd.read_csv(StringIO(session.get('gtex_df')))
        svs_df = pd.read_csv(StringIO(session.get('svs_df')))
        gwascat = pd.read_csv(StringIO(session.get('gwascat_df')))

        # Set up tables
        df_sorted = df.sort_values(by='Pvalue', ascending=False)
        df_sorted['Locus'] = df_sorted.apply(lambda x: str(x['Chrom']) + ':' + str(x['Position']), axis=1)
        df_sorted['Alleles'] = df_sorted.apply(lambda x: str(x['Ref']) + '/' + str(x['Alt']), axis=1)
        df_sorted = df_sorted[['Locus', 'Rsid', 'Gwas', 'Alleles', 'Pvalue']]
        df_sorted['Pvalue'] = df_sorted['Pvalue'].round(2)
        snps_table = df_sorted.to_dict(orient='records')
        svs_table = svs_df.to_dict(orient='records')
        gwas_table = gwascat.to_dict(orient='records')

        # Plot
        img = scatterplot(df = df, chrom = chrom, start_pos = start_pos, end_pos = end_pos, gwas = gwas, genes = genes, svs = svs, browse_type = browse_type, refGen = refGen, recomb_data = recomb_data, exons = exons, browse = browse, plotype = plotype)
        # set the plot url for showing on the application
        plot_url = base64.b64encode(img.getvalue()).decode()
        # then generate the GTEx plot
        img_gtex = gtex_heatmap(gtex_df)
        # set the plot url for showing on the application
        plot_gtex = base64.b64encode(img_gtex.getvalue()).decode()

        # return the html template and the url to the plot
        return render_template("exploration.html", plot_url=plot_url, plot_gtex=plot_gtex, table_snps=snps_table, table_svs=svs_table, table_gwas=gwas_table, browse_value=browse, gwas=gwas)
    else:
        gwas = []
        return render_template("exploration.html")

# Download SNP table
@app.route('/downloadSNPs')
def download_snp_table():
    if 'df' in session:
        # retrieve table data from session object
        csv_data = session['df']
        # convert CSV string to file-like object
        file_like = StringIO(csv_data)
        # create response with file-like object as download
        response = make_response(file_like.getvalue())
        response.headers["Content-Disposition"] = "attachment; filename=SNPs_table.csv"
        response.headers["Content-type"] = "text/csv"    
        return response
    return 'Data not available. Please go back and generate a plot before downloading.'

# Download SV table
@app.route('/downloadSVs')
def download_sv_table():
    if 'svs_df' in session:
        # retrieve table data from session object
        csv_data = session['svs_df']
        # convert CSV string to file-like object
        file_like = StringIO(csv_data)
        # create response with file-like object as download
        response = make_response(file_like.getvalue())
        response.headers["Content-Disposition"] = "attachment; filename=SVs_table.csv"
        response.headers["Content-type"] = "text/csv"    
        return response
    return 'Data not available. Please go back and generate a plot before downloading.'

# Download GWAS-Catalog table
@app.route('/downloadGWAScat')
def download_gwascat_table():
    if 'gwascat_df' in session:
        # retrieve table data from session object
        csv_data = session['gwascat_df']
        # convert CSV string to file-like object
        file_like = StringIO(csv_data)
        # create response with file-like object as download
        response = make_response(file_like.getvalue())
        response.headers["Content-Disposition"] = "attachment; filename=GWASCatalog_table.csv"
        response.headers["Content-type"] = "text/csv"    
        return response
    return 'Data not available. Please go back and generate a plot before downloading.'

# Download GTEx Expression
@app.route('/downloadGTExExpression')
def download_gtex_expression():
    if 'gtex_df' in session:
        # retrieve table data from session object
        csv_data = session['gtex_df']
        # convert CSV string to file-like object
        file_like = StringIO(csv_data)
        # create response with file-like object as download
        response = make_response(file_like.getvalue())
        response.headers["Content-Disposition"] = "attachment; filename=GTEx_expression.csv"
        response.headers["Content-type"] = "text/csv"    
        return response
    return 'Data not available. Please go back and generate a plot before downloading.'

# Download SNP plot
@app.route('/download_plot/SNP_plot.png')
def download_snp_plot():
    if 'df' in session:
        chrom = session['chrom']
        start_pos = session['start_pos']
        end_pos = session['end_pos']
        gwas = session['gwas']
        genes = session['genes']
        svs = session['svs']
        browse_type = session['browse_type']
        refGen = session['refGen']
        recomb_data = session['recomb_data']
        exons = session['exons']
        browse = session['browse']
        plotype = session['plotype']
        df = pd.read_csv(StringIO(session.get('df')))
        gtex_df = pd.read_csv(StringIO(session.get('gtex_df')))
        svs_df = pd.read_csv(StringIO(session.get('svs_df')))
        gwascat = pd.read_csv(StringIO(session.get('gwascat_df')))
        filename = 'SNP_plot.png'
        # Plot
        img = scatterplot(df = df, chrom = chrom, start_pos = start_pos, end_pos = end_pos, gwas = gwas, genes = genes, svs = svs, browse_type = browse_type, refGen = refGen, recomb_data = recomb_data, exons = exons, browse = browse, plotype = plotype)
        return send_file(img, mimetype='image/png', as_attachment=True, download_name=filename)
    return 'Data not available. Please go back and generate a plot before downloading.'

# Download SNP plot
@app.route('/download_gtex/GTEx_plot.png')
def download_gtex_plot():
    if 'df' in session:
        gtex_df = pd.read_csv(StringIO(session.get('gtex_df')))
        filename = 'GTEx_plot.png'
        # Plot
        img_gtex = gtex_heatmap(gtex_df)
        return send_file(img_gtex, mimetype='image/png', as_attachment=True, download_name=filename)
    return 'Data not available. Please go back and generate a plot before downloading.'

# About tab
@app.route('/about/')
def about():
    return render_template('about.html')

# User tab
@app.route('/user/')
def user():
    if 'user_id' in session:
        return redirect(url_for('dashboard'))
    if request.method == 'POST':
        if request.form['action'] == 'login':
            # Handle login form submission
            return 'Login button pressed'
        elif request.form['action'] == 'register':
            #Handle registration form submission
            return 'Register button pressed'
    else:
        # Render user.html template with GET request
        return render_template('user.html')

# Logout
@app.route('/logout')
def logout():
    # remove the username from the session if it is there
    session.pop('user_id', None)
    return redirect(url_for('login'))

# Annotation section
@app.route('/annotation/', methods=["GET", "POST"])
def annotation():
    #if 'user_id' in session:
    #    message = '%s You will find the results in your dashboard.' %(session['success_login'])
    #else:
    #    message = 'Please enter your email address below. The results will be sent to you:'
    #    email = '<textarea name="emailAnnot" style="width: 100%; height: 30px;" placeholder="Enter your email here"></textarea>'
    # check if the user has logged in -- in this case, show a messages, otherwise, show a textarea for the user to input the email address
    # email
    email = request.form.get('email', '')
    email = email.replace(' ', '').rstrip()
    # check if email is valid
    emailResponse = is_valid_email(email)
    # modify the message linked to submission
    messageSubmission = ''
    if request.method == "POST":
        # Read inputs
        # list of variants
        textarea_input = request.form.get('SNPlist', '')
        my_list = textarea_input.replace(' ', '').replace(';', '\n').replace(',', '\n').split('\n')
        # input type
        inpType = request.form["inputType"]
        # if the input type is not rsid, then read also the reference genome
        refGeno = request.form["annotRefGen"] if inpType != 'RsID' else 'GRCh37'
        # analysis type
        analType = request.form["analType"]
        # if the analysis type is GSEA, then look at the enrichment sources
        gsea_source = request.form.getlist("EnrichSource") if analType == 'GSEA' else []
        # QTL tissues
        qtl_tissues = request.form.getlist('QTLtissuesAnnot')
        # update message
        if emailResponse == True:
            messageSubmission = 'Your job has been submitted! Check your email shortly.'
            # run the script here
            command = run_annotation(my_list, inpType, refGeno, analType, gsea_source, qtl_tissues, email)
            print(command)
            command_run = subprocess.Popen(command, shell=True)
        else:
            messageSubmission = 'Email is not correct. Please check!'
    return render_template('annotation.html', validEmail=emailResponse, messageSubmission=messageSubmission) 

# Help tab
@app.route('/help/')
def help():
    return render_template('help.html')

# Download tab
@app.route('/download/', methods=["GET", "POST"])
def download():
    run_id = None
    messageError = None
    messageToUser = None
    if request.method == 'POST':
        run_id = request.form["run_id"]
        messageError, messageToUser = check_runID(run_id)
    return render_template('download.html', run_id=run_id, messageError=messageError, messageToUser=messageToUser)

# Actual download button for the annotation results
@app.route('/downloadResults/<run_id>', methods=["GET"])
def download_results(run_id):
    # check here the run ID
    filename = f"AnnotateMe_results_{run_id}.tar.gz"
    filepath = os.path.join("/Annotation/RUNS", filename)
    return send_file(filepath, as_attachment=True, download_name=filename)

# Global functions
# function to check whether the run_id is correct
def check_runID(run_id):
    run_id = run_id.replace(' ', '')
    if run_id == '':
        messageError = True
        messageToUser = 'Please insert a valid Run ID'
    else:
        # check whether the file exists
        if os.path.exists('/Annotation/RUNS/AnnotateMe_results_%s.tar.gz' %(run_id)):
            messageError = False
            messageToUser = 'Valid Run ID. Download will start soon'
        else:
            messageError = True
            messageToUser = 'Not valid Run ID. Try again.'
    return messageError, messageToUser

# function to extract coordinates based on input type
def readBrowseOption(data_path, browse, window, refGen):
    # check if input is in the form 1:100000
    if ':' in browse and '-' not in browse:
        chrom, start_pos, end_pos, browse_type = int(browse.split(':')[0]), int(browse.split(':')[1]) - window, int(browse.split(':')[1]) + window, 'Single position'
    elif ':' in browse and '-' in browse:
        chrom, start_pos, end_pos, browse_type = int(browse.split(':')[0]), int(browse.split(':')[1].split('-')[0]) - window, int(browse.split(':')[1].split('-')[1]) + window, 'Interval'
    elif 'rs' in browse:
        snps = pd.DataFrame(get_variants_by_variant_id(browse.replace(' ', '')).locations)
        chrom, pos, browse_type = snps['chromosomeName'][0], snps['chromosomePosition'][0], 'RsID'
        if refGen == 'GRCh37':
            converter = get_lifter('hg38', 'hg19')
            start_pos, end_pos = converter[chrom][pos][0][1] - window, converter[chrom][pos][0][1] + window
        else:
            start_pos, end_pos = pos - window, pos + window
    else:
        data_path = '%s/databases/Genes/genes_hg19.txt.gz' %(data_path) if refGen == 'GRCh37' else '%s/databases/Genes/genes_hg38.txt.gz' %(data_path)
        genes = list(os.popen('zgrep -i -w %s %s' %(browse.replace(' ', ''), data_path)))[0].split('\t')
        if refGen == "GRCh37":
            chrom, start_pos, end_pos, browse_type = genes[2].replace('chr', ''), int(genes[4]) - window, int(genes[5]) + window, 'Gene'
        else:
            chrom, start_pos, end_pos, browse_type = genes[1].replace('chr', ''), int(genes[3]) - window, int(genes[4]) + window, 'Gene'
    return chrom, start_pos, end_pos, browse_type

# function to run annotation analysis
def run_annotation(my_list, inpType, refGeno, analType, gsea_source, qtl_tissues, email):
    # Sample a number for randomization
    random_number = random.randint(1, 10000)
    # Take snps and save them
    fpath_server = '/Annotation/RUNS/annotateMe_input_%s.txt' %(random_number)
    filename = 'annotateMe_input_%s.txt' %(random_number)
    with open(fpath_server, 'w') as finp:
        for x in my_list:
            finp.write('%s\n' %(x))
    # Define log file
    log_filename_server = '/Annotation/RUNS/annotateMe_run_%s.log' %(random_number)
    # Take input type
    ftype = str(inpType)
    if ftype == 'Colon':
        ftype = 1
    elif ftype == 'Tab':
        ftype = 2
    else:
        ftype = 3
    # Take analysis type
    analysis_type = str(analType)
    analysis_type = 'mapping' if analysis_type == 'Annot' else 'enrichment'
    print(gsea_source)
    analysis_mode = ','.join(list(gsea_source)) if analysis_type == 'enrichment' else 'None'
    # Take gtex tissues
    gtex_tissues = ','.join(list(qtl_tissues))
    # Then run annotate me externally in background -- this depends on the analysis_type requested
    command = "Rscript /Annotation/BIN/MAIN.R %s %s %s %s %s %s %s %s > %s" %(filename, ftype, email, analysis_type, analysis_mode, gtex_tissues, refGeno, random_number, log_filename_server)
    return command

# function to extract GTEx information
def get_gtex(data_path, genes, refGen):
    genename_index = 0 if refGen == 'GRCh37' else -2
    if len(genes) >0:
        gtex_genes = []
        for x in genes:
            try:
                tmp = [x.rstrip().split('\t') for x in list(os.popen('zgrep %s %s/databases/GTEx_Analysis_2017_06_05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz' %(x[genename_index], data_path)))][0]
                gtex_genes.append(tmp)
            except:
                pass
    # convert to dataframe
    colnames = ["Name", "Description", "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Adrenal Gland", "Artery - Aorta", "Artery - Coronary", "Artery - Tibial", "Bladder", "Brain - Amygdala", "Brain - Anterior cingulate cortex (BA24)", "Brain - Caudate (basal ganglia)", "Brain - Cerebellar Hemisphere", "Brain - Cerebellum", "Brain - Cortex", "Brain - Frontal Cortex (BA9)", "Brain - Hippocampus", "Brain - Hypothalamus", "Brain - Nucleus accumbens (basal ganglia)", "Brain - Putamen (basal ganglia)", "Brain - Spinal cord (cervical c-1)", "Brain - Substantia nigra", "Breast - Mammary Tissue", "Cells - Cultured fibroblasts", "Cells - EBV-transformed lymphocytes", "Cervix - Ectocervix", "Cervix - Endocervix", "Colon - Sigmoid", "Colon - Transverse", "Esophagus - Gastroesophageal Junction", "Esophagus - Mucosa", "Esophagus - Muscularis", "Fallopian Tube", "Heart - Atrial Appendage", "Heart - Left Ventricle", "Kidney - Cortex", "Kidney - Medulla", "Liver", "Lung", "Minor Salivary Gland", "Muscle - Skeletal", "Nerve - Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)", "Small Intestine - Terminal Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole Blood"]
    gtex = pd.DataFrame(gtex_genes, columns=colnames)
    return gtex

# function to get data to plot given the gwas name, chromosome, start and end position
def get_data_plot(data_path, gwas, chrom, start_pos, end_pos, refGen):
    # set reference genome
    refPrefix = '' if refGen == 'GRCh37' else '_hg38'
    # iterate over the gwas to plot
    data_list = []
    for gw in gwas:
        # define command for tabix
        cmd = 'tabix %s/%s/chr%s_%s%s.txt.gz %s:%s-%s' %(data_path, gw, str(chrom), gw, refPrefix, str(chrom), str(start_pos), str(end_pos))
        # extract summary statistics and add gwas name (in case of the GWAS catalog, this is slightly different as the GWAS trait is included in the data)
        sumstats = [[gw] + x.rstrip().split('\t') for x in os.popen(cmd)] if gw != 'GWAS-Cat' else [x.rstrip().split('\t') for x in os.popen(cmd)]
        # append to data_list
        data_list = sumstats if len(data_list) == 0 else data_list + sumstats
    # convert data_list to dataframe -- make a slight difference for the GWAS-Catalog
    colnames = ['Gwas', 'Position', 'Chrom', 'Pvalue', 'Position_hg38', 'Rsid', 'Maf', 'Ref', 'Alt'] if refGen == 'GRCh37' else ['Gwas', 'Chrom', 'Pvalue', 'Position', 'Rsid', 'Maf', 'Ref', 'Alt']
    df = pd.DataFrame(data_list, columns=colnames)
    # make position and pvalue numeric
    df["Position"] = pd.to_numeric(df['Position'], errors='coerce')
    df["Pvalue"] = pd.to_numeric(df['Pvalue'], errors='coerce')
    # adjust the SNPs with pvalue 0
    zero_rows = df['Pvalue'] == 0
    df.loc[zero_rows, 'Pvalue'] = 5.00e-30
    # then take the -log10 of the pvalue
    df['Pvalue'] = -np.log10(df['Pvalue'])
    # then take the position value as megabases
    df["Position_plot"] = df["Position"] / 1000000
    # drop duplicates based on position and gwas
    df = df.drop_duplicates(subset=['Gwas', 'Position'])
    return df

# function to extract genes to plot
def extract_genes(data_path, chrom, start_pos, end_pos, refGen):
    # set prefix of the reference used
    refPrefix = 'hg19' if refGen == 'GRCh37' else 'hg38'
    # use tabix to find genes -- enlarge window by 50kb up and down
    cmd = 'tabix %s/databases/Genes/genes_%s.txt.gz chr%s:%s-%s' %(data_path, refPrefix, str(chrom), str(start_pos - 1000), str(end_pos + 1000))
    genes = [x.rstrip().split('\t') for x in os.popen(cmd)]
    # in case of GRCh38, there may be multiple entries for the same gene -- keep the one where the gene is largest
    if refGen == 'GRCh38':
        unique_genes = {}
        for g in genes:
            key = g[-1]
            if key not in unique_genes:
                unique_genes[key] = g
            else:
                if g[0] > unique_genes[key][0]:
                    unique_genes[key] = g
        genes = list(unique_genes.values())
    # exclude long non coding stuff
    if refGen == 'GRCh37':
        genes = [x for x in genes if 'LOC' not in x[0]]
    else:
        genes = [x for x in genes if 'LOC' not in x[-1]]
    # add y axis
    genes = [[*sublist, i] for i, sublist in enumerate(genes)]
    return genes

# function to extract genes to plot
def extract_sv(data_path, chrom, start_pos, end_pos, selected_sv_source, refGen):
    # set reference prefix
    refPrefix = 'hg19' if refGen == 'GRCh37' else 'hg38'
    # use tabix to find genes -- enlarge window by 50kb up and down
    cmd = 'tabix %s/databases/Structural_variants/str_set_%s_homogeneous.txt.gz chr%s:%s-%s' %(data_path, refPrefix, str(chrom), str(start_pos), str(end_pos))
    svs = [x.rstrip().split('\t') for x in os.popen(cmd)]
    # select based on the input selected
    svs = [x for x in svs if x[6] in selected_sv_source]
    # add y axis
    svs = [[*sublist, i] for i, sublist in enumerate(svs)]
    # convert to dataframe as well
    colnames = ['chrom', 'start', 'end', 'type', 'col', 'diff', 'source', 'y'] if refGen == 'GRCh37' else ['chrom', 'start', 'end', 'diff', 'type', 'col', 'source', 'y']
    svs_df = pd.DataFrame(svs, columns=colnames)
    return svs, svs_df

# function to extract recombination rates
def extract_recomb(data_path, chrom, start_pos, end_pos, refGen):
    # set reference prefix
    refPrefix = 'hg19' if refGen == 'GRCh37' else 'hg38'
    # use tabix to extract recombination rates
    cmd = 'tabix %s/databases/Recombination_rates/recombination_rates_%s.txt.gz chr%s:%s-%s' %(data_path, refPrefix, str(chrom), str(start_pos), str(end_pos))
    rec_rates = [x.rstrip().split('\t') for x in os.popen(cmd)]
    # convert to dataframe
    rec_rates_df = pd.DataFrame(rec_rates, columns=['chrom', 'position', 'rate', 'map'])
    rec_rates_df['position'] = pd.to_numeric(rec_rates_df['position'])/1000000
    rec_rates_df['rate'] = pd.to_numeric(rec_rates_df['rate'])
    return rec_rates_df

# function to extract GWAS catalog information
def extract_gwascatalog(data_path, chrom, start_pos, end_pos, refGen):
    # set reference prefix
    refPrefix = 'hg19' if refGen == 'GRCh37' else 'hg38'
    # use tabix to extract recombination rates
    cmd = 'tabix %s/databases/GWAS_catalog/Gwas_catalog_%s.txt.gz %s:%s-%s' %(data_path, refPrefix, str(chrom), str(start_pos), str(end_pos))
    gwascat = [x.rstrip().split('\t') for x in os.popen(cmd)]
    # convert to dataframe
    colnames = ['chr', 'pos', 'pos_hg38', 'snp', 'gene', 'trait', 'pubmed_id', 'pvalue'] if refGen == 'GRCh37' else ['chr', 'pos_hg37', 'pos', 'snp', 'gene', 'trait', 'pubmed_id', 'pvalue']
    gwascat = pd.DataFrame(gwascat, columns=colnames)
    # prepare table for output: add locus, remove redundant columns
    gwascat['Locus'] = gwascat.apply(lambda x: str(x['chr']) + ':' + str(x['pos']), axis=1)
    gwascat = gwascat[['Locus', 'snp', 'gene', 'trait', 'pubmed_id', 'pvalue']]
    gwas_table = gwascat.to_dict(orient='records')
    return gwas_table, gwascat

# function to draw the scatterplot given a dataframe and positions to plot
def scatterplot(df, chrom, start_pos, end_pos, gwas, genes, svs, browse_type, refGen, recomb_data, exons, browse, plotype):
    # define structure for the plot
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 12), gridspec_kw={'height_ratios': [2, 1.5, 1.5]})
    ##################
    # Scatterplot here
    # define colors and map them to datasets
    #all_colors = ['red', 'blue', 'orange', 'yellow', 'grey', 'green']
    all_colors = [x for x in list(webcolors.CSS3_NAMES_TO_HEX.keys()) if x != 'white']
    random.shuffle(all_colors)
    gwas_unique = df['Gwas'].unique().tolist()
    # sometimes using the GWAS-Cat data, there are too many colors needed and it crashes. In these cases, take the top 100 traits based on significance
    if len(gwas_unique) > 70:
        df_gwascat = df[df['Gwas'].isin(gwas_unique)]
        df_gwascat = df_gwascat[~df_gwascat['Gwas'].isin(gwas)]
        df_NOgwascat = df[df['Gwas'].isin(gwas)]
        df_gwascat = df[df['Gwas'].isin(gwas_unique)]
        # sort the dataframe by the numeric column in descending order
        df_sorted = df_gwascat.sort_values(by='Pvalue', ascending=False)
        # select the top 100 unique values from the character column
        unique_values = df_sorted['Gwas'].unique()[:70]
        # subset the dataframe based on the selected values
        subset = df_sorted[df_sorted['Gwas'].isin(unique_values)]
        # re-join data
        df = pd.concat([df_NOgwascat, subset])
        gwas_unique = df['Gwas'].unique().tolist()
    color_map = {val: all_colors[index] for index, val in enumerate(gwas_unique)}
    # Set subplot background color to grey
    ax1.set_facecolor('#f0f0f0')
    ax2.set_facecolor('#f0f0f0')
    ax3.set_facecolor('#f0f0f0')
    # Then put a white vertical grid
    minimum_value_pos = df['Position_plot'].min()
    maximum_value_pos = df['Position_plot'].max()
    num_values = 10     # define grid size
    equally_spaced_values = np.linspace(minimum_value_pos, maximum_value_pos, num_values)
    for val in equally_spaced_values:
        ax1.axvline(x=val, color='white', linestyle='-', linewidth=1.5, zorder = 1)
        ax2.axvline(x=val, color='white', linestyle='-', linewidth=1.5, zorder = 1)
        ax3.axvline(x=val, color='white', linestyle='-', linewidth=1.5, zorder = 1)
    # plot data depending on plot type
    if plotype == 'Scatter':
        colors = df['Gwas'].map(color_map)
        # scatter plot
        ax1.scatter(df['Position_plot'], df['Pvalue'], c=colors, zorder = 2)
    else:
        # smooth plot - need to iterate over gwas types
        for gw in gwas:
            if gw == 'GWAS-Cat':
                # if the whole GWAS-Catalog was selected, still use scatter as it is too sparse to do densities
                tmp = df.loc[~df['Gwas'].isin(gwas)].copy()
                colors = tmp['Gwas'].map(color_map)
                ax1.scatter(tmp['Position_plot'], tmp['Pvalue'], c=colors, zorder = 2)
            else:
                tmp = df.loc[df['Gwas'] == gw].copy()
                # use rolling function to select local maxima
                tmp['smoothed'] = tmp['Pvalue'].rolling(window=20).max()
                # fit lowess
                lowess = sm.nonparametric.lowess(tmp['smoothed'], tmp['Position_plot'], frac=0.02)
                # fill the area under the lowess line
                ax1.fill_between(lowess[:, 0], lowess[:, 1], color=color_map[gw], alpha=0.4)
                ax1.plot(lowess[:, 0], lowess[:, 1], color=color_map[gw], linewidth=3)
    # set y-lim -- max pvalue + 10% of the max value
    max_p_assoc = int(max(df['Pvalue'])) + int(max(df['Pvalue']))*0.08
    ax1.set_ylim(0, max(max_p_assoc, 6))
    # if rsid was requested, the corresponding marker will be in a different color (black)
    if browse_type == 'RsID':
        subset_rsid = df.loc[df['Rsid'] == browse.lower()].copy()
        # add colors
        colors_rsid = subset_rsid['Gwas'].map(color_map)
        # also label snp of interest
        subset_rsid['Position_label_x'] = subset_rsid['Position_plot'] - (end_pos - start_pos)/1000000*0.25
        subset_rsid['Position_label_y'] = subset_rsid['Pvalue']
        subset_rsid['Labels'] = 'Input SNP ~ ' + subset_rsid['Gwas']
        ax1.scatter(subset_rsid['Position_plot'], subset_rsid['Pvalue'], c=colors_rsid, marker = '^')
        subset_rsid_lst = subset_rsid.values.tolist()
        if len(gwas_unique) <5:
            for x in subset_rsid_lst:
                ax1.annotate(x[-1], xy=(x[-4], x[-2]), xytext=(x[-3], x[-2]), arrowprops=dict(facecolor='black', shrink=0.002, width=1, mutation_scale=10))
    # add text for rsid of the most significant variant: the selection is different in case of the GWAS-Cat: in this case we take only the top 3 traits
    top_snp = df.loc[df.groupby('Gwas')['Pvalue'].idxmax(), ['Position_plot', 'Pvalue', 'Gwas', 'Rsid']]
    if 'GWAS-Cat' in gwas:
        top_snp_noCat = top_snp[top_snp['Gwas'].isin(gwas)]
        top_snp_Cat = top_snp.sort_values('Pvalue', ascending=False)
        top_snp_Cat = top_snp_Cat.head(3)
        top_snp = pd.concat([top_snp_Cat, top_snp_noCat], axis=0)
    top_snp['Position_label_x'] = top_snp['Position_plot'] + (end_pos - start_pos)/1000000*0.04; top_snp['Position_label_y'] = top_snp['Pvalue']
    top_snp['Labels'] = top_snp['Rsid'] + ' ~ ' + top_snp['Gwas']; top_snp_lst = top_snp.values.tolist()
    for x in top_snp_lst:
        if browse_type != 'Rsid':
            ax1.annotate(x[-1], xy=(x[0], x[1]), xytext=(x[4], x[1]), arrowprops=dict(facecolor='black', shrink=0.002, width=1, mutation_scale=10))
        else:
            if x[3] not in str(subset_rsid['Rsid']):
                ax1.annotate(x[-1], xy=(x[0], x[1]), xytext=(x[4], x[1]), arrowprops=dict(facecolor='black', shrink=0.002, width=1, mutation_scale=10))
    # if requested, plot recombination rates
    if isinstance(recomb_data, pd.DataFrame):
        ax1_2 = ax1.twinx()
        ax1_2.plot(recomb_data['position'], recomb_data['rate'], c='green')
        ax1_2.set_ylabel('Recombination Rates (cM/Mb)', color='green')
        ax1_2.set_ylim(-0.5, 100)
    # add title
    title = '%s:%s-%s Input:%s (%s)' %(str(chrom), str(start_pos), str(end_pos), browse_type, refGen)
    ax1.set_title(title)
    # axes names
    ax1.set_xlabel("")
    ax1.set_ylabel("-log10(P-value)")
    # disable x-axis values
    ax1.set_xticks([])
    # Add the legend -- if the number of characters in the legend name is >20, then put a . and truncate it
    for key in list(color_map.keys()):
        if len(key) > 20:
            new_key = key[:20] + "."
            color_map[new_key] = color_map.pop(key)
    handles = [plt.plot([],[],color=color_map[color], marker="o", ls="", label=color)[0] for color in color_map.keys()]
    # plot the legend outside the plot (on top). Use columns containing 7 values. If more that 28 elements need to be plotted, make it extra small
    if len(color_map.keys()) >30:
        legend = ax1.legend(handles=handles, loc = 'upper center', bbox_to_anchor=(0.5, 1.45), fontsize=5.4, ncol = math.ceil(len(color_map.keys())/10))
    elif len(color_map.keys()) >15:
        legend = ax1.legend(handles=handles, loc = 'upper center', bbox_to_anchor=(0.5, 1.45), fontsize=7, ncol = math.ceil(len(color_map.keys())/8))
    elif len(color_map.keys()) >=10:
        legend = ax1.legend(handles=handles, loc = 'upper center', bbox_to_anchor=(0.5, 1.45), fontsize=8.5, ncol = math.ceil(len(color_map.keys())/8))
    else:
        legend = ax1.legend(handles=handles, loc = 'upper center', bbox_to_anchor=(0.5, 1.35), fontsize=10, ncol = math.ceil(len(color_map.keys())/8))
    ##################
    # Gene track here
    max_y = max(g[-1] for g in genes) if len(genes) >0 else 1
    index_pos = 4 if refGen == 'GRCh37' else 5
    index_text = 0 if refGen == 'GRCh37' else -2
    index_exon = 9 if refGen == 'GRCh37' else 8
    if len(genes) >0:
        for g in genes:
            ax2.plot([int(g[index_pos])/1000000, int(g[index_pos + 1])/1000000], [int(g[-1]), int(g[-1])], 'r-', linewidth=3)
            if len(genes) < 10:
                ax2.text(int(g[index_pos])/1000000 + (int(g[index_pos + 1])/1000000 - int(g[index_pos])/1000000)/2, int(g[-1])+max_y*0.06, g[index_text], ha='center', fontsize = 10)
            elif len(genes) <20:
                ax2.text(int(g[index_pos])/1000000 + (int(g[index_pos + 1])/1000000 - int(g[index_pos])/1000000)/2, int(g[-1])+max_y*0.03, g[index_text], ha='center', fontsize = 7)
            # show exons if requested
            if exons == 'Yes':
                # extract exon starts and end
                exons_start = [int(x)/1000000 for x in g[index_exon].split(',') if x != '']
                exons_end = [int(x)/1000000 for x in g[index_exon+1].split(',') if x != '']
                ax2.plot([exons_start, exons_end], [int(g[-1]), int(g[-1])], 'r-', linewidth=5)
    else:
        ax2.text((start_pos + (end_pos-start_pos)/2)/1000000, 0.5, "No Genes to plot", ha='center')
    # axes names
    ax2.set_xlabel("")
    ax2.set_ylabel("Gene track")
    # disable axes values
    ax2.set_xticks([])
    ax2.set_yticks([])
    # set x-limit
    ax2.set_xlim(int(start_pos)/1000000, int(end_pos)/1000000)
    ax2.set_ylim(0-max_y*0.20, max_y+max_y*0.20)
    ##################
    # Structural variants
    max_y = max(sv[-1] for sv in svs) if len(svs) >0 else 1
    index_col = 4 if refGen == 'GRCh37' else 5
    index_sv = 3 if refGen == 'GRCh37' else 4
    index_diff = -3 if refGen == 'GRCh37' else -5
    if len(svs) >0:
        for sv in svs:
            # if the size is too small (<5% of the window), then plot a square
            if int(sv[index_diff]) < (end_pos - start_pos)*0.001:
                ax3.plot(int(sv[1])/1000000, sv[-1], color=sv[index_col], marker = 's', markersize = 5)
            else:
                ax3.plot([int(sv[1])/1000000, int(sv[2])/1000000], [sv[-1], sv[-1]], color=sv[index_col], linewidth=5)
    else:
        ax3.text((start_pos + (end_pos-start_pos)/2)/1000000, 0.5, "No Structural variants to plot", ha='center')
    # axes names
    ax3.set_xlabel("Genomic Position (Mb)")
    ax3.set_ylabel("Structural variants")
    # disable axes values
    ax3.set_yticks([])
    # set x-limit
    ax3.set_xlim(int(start_pos)/1000000, int(end_pos)/1000000)
    ax3.set_ylim(0-max_y*0.10, max_y+max_y*0.20)
    # Add the legend
    color_map_sv = {x[index_sv]: x[index_col] for x in svs}
    handles_sv = [plt.plot([],[],color=color_map_sv[color], marker="s", ls="", label=color)[0] for color in color_map_sv.keys()]
    if len(svs) >0:
        ax3.legend(handles=handles_sv, loc = 'upper center', ncol=6, fontsize=8, markerscale=0.75)
    #plt.show()         # enable this for debuggin
    # save image
    img = io.BytesIO()
    fig.savefig(img, format="png")
    img.seek(0)
    return img

# function to generate the GTEx heatmap
def gtex_heatmap(gtex_info):
    # select the columns to include in the heatmap
    cols_heatmap = gtex_info.columns[2:]
    tmp = gtex_info[cols_heatmap]
    labs_genes = list(gtex_info['Description'])
    labs_tissues = list(cols_heatmap)
    # convert all columns from the third on to numeric
    tmp = tmp.apply(pd.to_numeric)
    # plot the heatmap -- the size of the heatmap depends on the number of genes to be plotted
    if len(tmp) <= 5:
        fig, ax = plt.subplots(figsize=(16, 4))
    elif len(tmp) <= 10:
        fig, ax = plt.subplots(figsize=(16, 5))
    elif len(tmp) <= 15:
        fig, ax = plt.subplots(figsize=(16, 6))
    else:
        fig, ax = plt.subplots(figsize=(16, 7))
    im = ax.imshow(tmp, cmap='viridis', interpolation='nearest')
    # Set the x and y labels
    ax.set_xticks(np.arange(len(tmp.columns)))
    ax.set_yticks(np.arange(len(tmp.index)))
    ax.set_xticklabels(labs_tissues, fontsize=7)
    ax.set_yticklabels(labs_genes, fontsize=8)
    # Rotate the x labels if they overlap
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    # Set the colorbar
    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.5)
    # Set the title
    ax.set_title('GTEx expression')
    img_gtex = io.BytesIO()
    fig.savefig(img_gtex, format="png")
    img_gtex.seek(0)
    return img_gtex

# function to check if an email is valid
def is_valid_email(email):
    pattern = r'^[\w\.-]+@[\w\.-]+\.\w+$'
    return re.match(pattern, email) is not None

if __name__ == "__main__":
    app.run()