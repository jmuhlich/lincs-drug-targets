import sys
import os
import requests
import zipfile
import hashlib

url = 'http://www.drugbank.ca/system/downloads/current/drugbank.xml.zip'
datafile_name = 'drugbank.xml'
datafile_sha1 = '910619452ca3fe9659c9b893f368bca35ff268ab'

def verify_hash():
    try:
        datafile = open(datafile_name)
        sha1 = hashlib.sha1()
        for chunk in iter(lambda: datafile.read(2**13), b''): 
            sha1.update(chunk)
        if sha1.hexdigest() == datafile_sha1:
            return True
        else:
            print "Data file corrupt!"
    except IOError as e:
        print "Error opening data file!\n--> ", e
    return False

def fetch_datafile():
    tmpfile = os.tmpfile()
    req = requests.get(url, stream=True)
    tmpfile.writelines(req)
    req.close()
    tmpfile.flush()
    tmpfile.seek(0)
    zipdata = zipfile.ZipFile(tmpfile, 'r')
    zipdata.extract(datafile_name, path=os.path.dirname(sys.argv[0]))


if verify_hash():
    print "Data file is OK."
else:
    print "Downloading file..."
    fetch_datafile()
    print "Done."
