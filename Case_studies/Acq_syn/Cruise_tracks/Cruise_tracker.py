'''
This short script will download data from the EUMETSAT Data Store to match a 
supplied cruise track. All user parameters are passed to the script via a
configuration file (config.ini, by default). If you have not used this script
before, we recommend you work through the companion Jupyter Notebook
Cruise_tracker.ipynb for guidance.

License: MIT
Author: Ben Loveday
Copyright: EUMETSAT, 2024
'''

from argparse import ArgumentParser
import cartopy
import configparser
import datetime
import eumdac
import getpass
import numpy as np
import os
import pandas as pd
from pathlib import Path
import shutil
import warnings
import xarray as xr
import zipfile
import sys

warnings.filterwarnings('ignore')

parser = ArgumentParser()
parser.add_argument("-c", "--config_file", help="configuration file",
    default="config.ini")
parser.add_argument("--search_only", help="switch to search only (no download)",
    action='store_true')
args = parser.parse_args()

if __name__ == '__main__':

    # read configuration parameters
    config_file = args.config_file
    config = configparser.ConfigParser()
    config.read(config_file)

    # pre-process some variables
    tvar = config["cruise_processing"]["time_variable"]
    latvar = config["cruise_processing"]["latitude_variable"]
    lonvar = config["cruise_processing"]["longitude_variable"]
    spatial_tolerance = float(config["data_selection"]["spatial_tolerance"])
    temporal_tolerance = float(config["data_selection"]["temporal_tolerance"])

    if config["data_selection"]["download_data"] == "True":
        download_data = True
    else:
        download_data = False

    components = config["data_selection"]["components"].split(',')
    flags = config["data_flagging"]["flags"].split(',')

    # Create a download directory for our products
    download_dir = os.path.join(os.getcwd(), config["data_selection"]["product_dir"])
    os.makedirs(download_dir, exist_ok=True)

    # process cruise data
    cruise_track_file = config["cruise_processing"]["cruise_file"]
    df = pd.read_csv(cruise_track_file,
        header=int(config["cruise_processing"]["header_lines"]),
        delimiter=config["cruise_processing"]["delimiter"])

    df_subset = df.copy()
    if int(config["cruise_processing"]["max_index"]) != -1:
        df_subset = df_subset[0:int(config["cruise_processing"]["max_index"])]
    if int(config["cruise_processing"]["min_index"]) != -1:
        df_subset = df_subset[int(config["cruise_processing"]["min_index"]):]
    if int(config["cruise_processing"]["stride"]) != -1:
        df_subset = df_subset[0::int(config["cruise_processing"]["stride"])]    

    # load eumdac credentials
    eumdac_credentials_file = Path(Path.home() / '.eumdac' / 'credentials')

    if os.path.exists(eumdac_credentials_file):
        consumer_key, consumer_secret = Path(eumdac_credentials_file).read_text().split(',')
    else:
        # creating authentication file
        consumer_key = input('Enter your consumer key: ')
        consumer_secret = getpass.getpass('Enter your consumer secret: ')
        try:
            os.makedirs(os.path.dirname(eumdac_credentials_file), exist_ok=True)
            with open(eumdac_credentials_file, "w") as f:
                f.write(f'{consumer_key},{consumer_secret}')
        except:
            pass
            
    token = eumdac.AccessToken((consumer_key, consumer_secret))
    print(f"This token '{token}' expires {token.expiration}")

    # instantiate data store
    datastore = eumdac.DataStore(token)

    # search for matching data
    collection = datastore.get_collection(config["data_selection"]["collectionID"])

    product_list = []
    product_index = []
    product_times = []
    product_lats = []
    product_lons = []

    for index, row in df_subset.iterrows():
        print(f"[{str(index).zfill(5)}] Checking available products for:\t{row[tvar]}\t\tLat: {row[latvar]},\t\tLon: {row[lonvar]}")

        # make a very small box around our cruise point
        ROI = [[row[lonvar] - spatial_tolerance, row[latvar] - spatial_tolerance],
               [row[lonvar] - spatial_tolerance, row[latvar] + spatial_tolerance],
               [row[lonvar] + spatial_tolerance, row[latvar] + spatial_tolerance],
               [row[lonvar] + spatial_tolerance, row[latvar] - spatial_tolerance],
               [row[lonvar] - spatial_tolerance, row[latvar] - spatial_tolerance]]

        # convert this to a WKT polygon
        polygon = 'POLYGON(({}))'.format(','.join(["{} {}".format(*coord) for coord in ROI]))

        # make a very small time window around our cruise point time
        dtstart = datetime.datetime.strptime(row[tvar], config["cruise_processing"]["tformat"]) - datetime.timedelta(seconds=temporal_tolerance/2)
        dtend = datetime.datetime.strptime(row[tvar], config["cruise_processing"]["tformat"]) + datetime.timedelta(seconds=temporal_tolerance/2)

        # find the matching products
        products = collection.search(dtstart=dtstart, dtend=dtend, geo=polygon, timeliness=config["data_selection"]["timeliness"])
        for product in products:
            print(f"        Found: {product}")
            product_list.append(product)
            product_index.append(index)
            product_times.append(row[tvar])
            product_lats.append(row[latvar])
            product_lons.append(row[lonvar])

    print(f"Identified {len(np.unique(product_list))} unique products for download")

    # write out matching files
    output_file = os.path.join(download_dir, config["cruise_processing"]["cruise_file"].split('.')[0] + "_output.csv")

    if os.path.exists(output_file):
        os.remove(output_file)

    with open(output_file, 'a') as f:
        f.write("index,time,latitude,longitude,product\n")
        for index, product_time, product_lat, product_lon, product in zip(product_index, product_times, product_lats, product_lons, product_list):
            f.write(f'{index},{product_time},{product_lat},{product_lon},{product}\n')

    # backup cruise and config files
    shutil.copyfile(cruise_track_file, os.path.join(download_dir, os.path.basename(cruise_track_file)))
    shutil.copyfile(config_file, os.path.join(download_dir, os.path.basename(config_file)))

    if args.search_only:
        print("Searching only....complete")
        sys.exit()

    # downloading matching files
    if download_data:
        for product, count in zip(np.unique(product_list), range(len(np.unique(product_list)))):
        
            if "all" in components:
        
                # download the required products
                with product.open() as fsrc, open(os.path.join(download_dir, fsrc.name), mode='wb') as fdst:
                    print(f'Downloading ({count+1}/{len(np.unique(product_list))}) {fsrc.name}.')
                    shutil.copyfileobj(fsrc, fdst)
            
                # Unzip the required products
                with zipfile.ZipFile(fdst.name, 'r') as zip_ref:
                    for file in zip_ref.namelist():
                        if file.startswith(str(product)):
                            zip_ref.extract(file, download_dir)
            
                # tidy up the required products
                os.remove(fdst.name)
        
            else:
                
                product_download_directory = os.path.join(download_dir, str(product))
                os.makedirs(product_download_directory, exist_ok=True)
        
                # download the required product components
                for entry in product.entries:
                    res = [ele for ele in components if(ele == os.path.basename(entry))]
                    if res:
                        with product.open(entry=entry) as fsrc, open(os.path.join(product_download_directory, fsrc.name),
                                                                    mode='wb') as fdst:
                            print(f'Downloading ({count+1}/{len(np.unique(product_list))}) {product}: {fsrc.name}.')
                            shutil.copyfileobj(fsrc, fdst)

    print('Complete')