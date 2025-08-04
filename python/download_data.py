#############################
### Download data imports ###

import os
import argparse
import requests
import tarfile
import gzip
import zipfile




# General function to download files from a URL and optionally extract them
def download_file_url(url: str, 
                      outdir: str, 
                      extract_gz: bool=False, 
                      overwrite: bool=False, 
                      extract_zip: bool=False, 
                      extract_tar: bool=False):
    """
    Will download file from url to outdir/filename
    filename is generated from the last portion of the url split by "/".
    Extract of compressed data can optionally be performed for tar, gz, and zip
    """
    
    # Format filepath and name
    filepath = os.path.join(outdir, url.split("/")[-1])
    fname = filepath.split("/")[-1]
    
    # Check if file already exists
    if os.path.isfile(filepath) and (overwrite == False):
        print("- Warning, file {} already exists... Set overwrite to True to download and replace".format(filepath))
        return
    
    
    # Download file and write in chuncks
    r = requests.get(url, stream=True)
    with open(filepath, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
        exit_statement = ["- Download of {} succesful...".format(url)]
    
    # Extract tar
    if extract_tar != False:
        file = tarfile.open(filepath)
        file.extractall(outdir)
        file.close()
        exit_statement.append("- Tar extract of {} succesful...".format(filepath))


    # Extract gzip
    elif extract_gz != False:
        outpath = os.path.join(outdir, fname.replace(".gz", ""))
        buffer_size = 1024 * 1024
        with gzip.open(filepath, 'rb') as f_in, open(outpath, 'wb') as f_out:
            while True:
                chunk = f_in.read(buffer_size)
                if not chunk:
                    break
                f_out.write(chunk)
        exit_statement.append("- Gzipped extract of {} successful...".format(filepath))
    
    elif extract_zip != False:
        with zipfile.ZipFile(filepath, 'r') as zip_ref:
            zip_ref.extractall(outdir)
        exit_statement.append("- Zipped extract of {} successful...".format(filepath))
    
    print("\n".join(exit_statement))


if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Downloads monarch data necessary for running pheval-semphen')
        parser.add_argument("-d", "--data_dir", help="Directory to write files to", required=True, type=str)
        return parser.parse_args()

    args = parse_input_command()
    ############################

    # Download directory
    download_dir = args.data_dir

    # Download links
    phenio_url = 'http://data.monarchinitiative.org/monarch-kg/latest/phenio.db.gz'
    monarch_kg_url = 'https://data.monarchinitiative.org/monarch-kg-dev/latest/monarch-kg.tar.gz'
    hgnc_mapping_url = 'https://data.monarchinitiative.org/mappings/latest/gene_mappings.sssom.tsv'

    # Create semphen data directory if does not exist
    if not os.path.isdir(download_dir):
        print("- Creating project directory at {}".format(download_dir))
        os.makedirs(download_dir, exist_ok=True)

    # Download data
    print("- Downloading necessary data files to {}".format(download_dir))
    download_file_url(phenio_url, download_dir, extract_gz=True, extract_zip=False, extract_tar=False, overwrite=False)
    download_file_url(monarch_kg_url, download_dir, extract_gz=False, extract_zip=False, extract_tar=True, overwrite=False)
    download_file_url(hgnc_mapping_url, download_dir, extract_gz=False, extract_zip=False, extract_tar=False, overwrite=False)