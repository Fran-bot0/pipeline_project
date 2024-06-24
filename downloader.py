#!/usr/bin/env python

from bs4 import BeautifulSoup
import os
import sys
import requests
import urllib.error
import urllib.request


BASE_URL = 'https://ftp.ensembl.org/pub/release-112/'

organism = sys.argv[1]

def get_href(url, name=organism):
    '''Gets the necessary href for the download'''
    response =  requests.get(url).text
    soup = BeautifulSoup(response, 'html.parser')
    for a in soup.find_all('a', href=True):
        if 'dna.primary_assembly.fa.gz' in a['href']:
            href = a['href']
            return href
        elif 'dna.toplevel.fa.gz' in a['href']:
            href = a['href']
            return href
        elif ('112.gtf.gz' in a['href']) or (f'{name}.vcf.gz' in a['href']):
            href = a['href']
            return href
    return None


def file_download(d_url, file_name):
    '''Takes an url and a file name. Downloads the file and names it.'''
    if d_url:
        if not os.path.isfile(file_name):
            try:
                urllib.request.urlretrieve(d_url, file_name)
            except urllib.error.HTTPError as e:
                print(f'{e}. It was not possible to download the file {file_name}.')
        else:
            print(f'{file_name} is already in this directory. Skipping this download')
    else:
        print(f'The {file_name} file you are trying to download is not available in Ensembl FTP Download')


primary_assembly_url = f'{BASE_URL}fasta/{organism}/dna/'
vcf_url = f'{BASE_URL}variation/vcf/{organism}/'
gtf_url = f'{BASE_URL}gtf/{organism}/'


primary_assembly_href = get_href(url=primary_assembly_url)
vcf_href = get_href(url=vcf_url)
gtf_href = get_href(url=gtf_url)


file_download(d_url=primary_assembly_url+primary_assembly_href, file_name=f'primary_assembly_{organism}.fa.gz')
file_download(d_url=vcf_url+vcf_href, file_name=f'vcf_{organism}.vcf.gz')
file_download(d_url=gtf_url+gtf_href, file_name=f'gtf_{organism}.gtf.gz')
