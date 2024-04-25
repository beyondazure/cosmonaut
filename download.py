# Written by beyondazure (Daniel Kabanovsky) at Boston College

import base64
import requests
import sys
import os

def download(file, description, version, auth_string):
    # encoding user credentials as base64
    string_bytes = auth_string.encode('ascii')
    base64_bytes = base64.b64encode(string_bytes)
    authorization_string = base64_bytes.decode('ascii')

    url = "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted"

    # setting up parameters
    params = {
        "path": version + "/cosmic/v99/" + file,
        "bucket": "downloads"
    }
    headers = {
        "Authorization": "Basic " + authorization_string
    }

    # making the GET request
    response = requests.get(url, params=params, headers=headers)

    if response.status_code == 200:
        # continue
        data = response.json()
        new_url = data.get("url")
        output = description
        # creating the output file
        if os.path.exists(output):
            postfix = 1
            while os.path.exists(output):
                output = output + "_" + str(postfix) + ".tar"
                postfix += 1
        else:
            output = description + ".tar"

        # formulating the curl request
        curl_request = "curl \"" + new_url + "\" --output " + output

        # executing the curl request and expanding / unzipping the downloaded file if it succeeded
        if os.system(curl_request) == 0:
            print("Success: downloaded " + file)
            tar_cmd = "tar -xf " + output
            if os.system(tar_cmd) != 0:
                print("Error while expanding file.")
                sys.exit()
            else:
                os.system("rm -f " + output)
                filename = file.split(".")[0] + ".tsv.gz"
                filename = filename.replace("Tsv_", "")
                unzip_cmd = "gunzip " + filename
                if os.system(unzip_cmd) != 0:
                    print("Error while unzipping file.")
                    sys.exit()
                else:
                    os.system("rm -f " + filename)
                    return filename.replace(".gz", "")
        # error handling
        else:
            print("Error while downloading file.")
            sys.exit()
    else:
        print("Error while downloading file: HTTP status code " + response.status_code)
        sys.exit()
    
def scripted_download():
    # obtaining/validating user input
    genomic_file = input("Copy the FULL name of the genomic screen TSV file here: ")
    targeted_file = input("Copy the FULL name of the targeted screen TSV file here: ")
    classification_file = input("Copy the FULL name of the classification TSV file here: ")
    version = input("Enter reference genome (GRCh37 / GRCh38) or leave blank (default GRCh38): ")
    if version not in ["grch37", "grch38"]:
        version = "grch38"
    else:
        version = version.lower()
    auth_string = input("Enter your authorization string: email@example.com:cosmic_password (this information is NOT stored by the program): ")
    files = [genomic_file, targeted_file, classification_file]
    descriptions = ["genomic", "targeted", "classification"]
    output = []
    # iterating over target files and invoking the download() method
    for file in files:
        new_file = download(file, descriptions[files.index(file)], version, auth_string)
        output.append(new_file)

    return output
