# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 17:46:25 2022

@author: Wenhao
"""


import time
import re
import json
import requests
import pandas as pd


headers = {
    'accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
    'accept-encoding': 'gzip, deflate, br, zstd',
    'accept-language': 'en-US,en;q=0.5',
    'Connection': 'keep-alive',
    'cookie': 'cf_clearance=YUnuM2C2Zzm9aGr2NiO0np9UA_04K.Djjyr6QsmpcaU-1762150990-1.2.1.1-hj0Fm3pLI7O986EqmMKCfbbC4ZgOE0XvtCEtz3hGcD_bW11dMBB7Yz8SYHHVoIXUiCM2yPedUBYfPBGmhKHBQ61HiVkqFX6Zpnvt2fP1to8oVz1vo.5fYqBITbVG8GdyefCHprKeBbvirZETu5vePzWkx.fLDhluV7cbGQc3vGVbFQOS6GS7JgLSOpXD.JtEqsAhGB6VU5_0uU6Xqzund5yRvv.AqSF.YP_7ScPaZlM',
    'Host': 'go.drugbank.com',
    'Priority': 'u=0, i',
    'Referer': 'https://go.drugbank.com/drugs/DB00317/clinical_trials/aggregate.json?&columns%5B0%5D%5Bdata%5D=0&columns%5B0%5D%5Bname%5D=&columns%5B0%5D%5Bsearchable%5D=true&columns%5B0%5D%5Borderable%5D=true&columns%5B0%5D%5Bsearch%5D%5Bvalue%5D=&columns%5B0%5D%5Bsearch%5D%5Bregex%5D=false&columns%5B1%5D%5Bdata%5D=1&columns%5B1%5D%5Bname%5D=&columns%5B1%5D%5Bsearchable%5D=true&columns%5B1%5D%5Borderable%5D=true&columns%5B1%5D%5Bsearch%5D%5Bvalue%5D=&columns%5B1%5D%5Bsearch%5D%5Bregex%5D=false&columns%5B2%5D%5Bdata%5D=2&c…5D=true&columns%5B3%5D%5Borderable%5D=true&columns%5B3%5D%5Bsearch%5D%5Bvalue%5D=&columns%5B3%5D%5Bsearch%5D%5Bregex%5D=false&columns%5B4%5D%5Bdata%5D=4&columns%5B4%5D%5Bname%5D=&columns%5B4%5D%5Bsearchable%5D=true&columns%5B4%5D%5Borderable%5D=true&columns%5B4%5D%5Bsearch%5D%5Bvalue%5D=&columns%5B4%5D%5Bsearch%5D%5Bregex%5D=false&start=0&length=50&search%5Bvalue%5D=&search%5Bregex%5D=false&__cf_chl_tk=chYY5G4WUDWNmN4.mn937cHZe.pMRcnxIvW00plsSBU-1762150985-1.0.1.1-kAMug8KwZ7p8GhYGG08zqDjGfSxTr8hRgm8MyPD_YNs',
    'Sec-Fetch-Dest': 'document',
    'Sec-Fetch-Mode': 'navigate',
    'Sec-Fetch-Site': 'same-origin',
    'Sec-GPC': '1',
    'Upgrade-Insecure-Requests': '1',
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:143.0) Gecko/20100101 Firefox/143.0'
}


def clinical_trial_of_drug(drug_id):
    error_url_ls = []
    start = 0
    length = 50
    df_ls = []
    while True:
        url = f'https://go.drugbank.com/drugs/{drug_id}/clinical_trials/aggregate.json?&columns%5B0%5D%5Bdata%5D=0&columns%5B0%5D%5Bname%5D=&columns%5B0%5D%5Bsearchable%5D=true&columns%5B0%5D%5Borderable%5D=true&columns%5B0%5D%5Bsearch%5D%5Bvalue%5D=&columns%5B0%5D%5Bsearch%5D%5Bregex%5D=false&columns%5B1%5D%5Bdata%5D=1&columns%5B1%5D%5Bname%5D=&columns%5B1%5D%5Bsearchable%5D=true&columns%5B1%5D%5Borderable%5D=true&columns%5B1%5D%5Bsearch%5D%5Bvalue%5D=&columns%5B1%5D%5Bsearch%5D%5Bregex%5D=false&columns%5B2%5D%5Bdata%5D=2&columns%5B2%5D%5Bname%5D=&columns%5B2%5D%5Bsearchable%5D=true&columns%5B2%5D%5Borderable%5D=true&columns%5B2%5D%5Bsearch%5D%5Bvalue%5D=&columns%5B2%5D%5Bsearch%5D%5Bregex%5D=false&columns%5B3%5D%5Bdata%5D=3&columns%5B3%5D%5Bname%5D=&columns%5B3%5D%5Bsearchable%5D=true&columns%5B3%5D%5Borderable%5D=true&columns%5B3%5D%5Bsearch%5D%5Bvalue%5D=&columns%5B3%5D%5Bsearch%5D%5Bregex%5D=false&columns%5B4%5D%5Bdata%5D=4&columns%5B4%5D%5Bname%5D=&columns%5B4%5D%5Bsearchable%5D=true&columns%5B4%5D%5Borderable%5D=true&columns%5B4%5D%5Bsearch%5D%5Bvalue%5D=&columns%5B4%5D%5Bsearch%5D%5Bregex%5D=false&start={start}&length={length}&search%5Bvalue%5D=&search%5Bregex%5D=false'
        
        response = requests.get(url, headers=headers)
        for i in range(10):
            if response.status_code == 200:
                break
            else:
                time.sleep(10)
                response = requests.get(url, headers=headers)
                
        if response.status_code != 200:
            print(f'Error! The request status code was {response.status_code}.')
            error_url_ls.append(url)
            continue
        
        data = json.loads(response.text) # dict_keys(['draw', 'recordsTotal', 'recordsFiltered', 'data'])
        PHASE = [re.search(r"'>([\w ]+?)</span>", e[0]).group(1) if '</span>' in e[0] else e[0] for e in data['data']]
        STATUS = [e[1] for e in data['data']]
        PURPOSE = [e[2] if '</span>' not in e[2] else re.search(r"'>([\w ]+?)</span>", e[2]).group(1) for e in data['data']]
        CONDITIONS = []
        CONDITIONS_LINK = []
        CONDITIONS_ls = [e[3].split(' / ') for e in data['data']]
        for e in CONDITIONS_ls:
            CONDITIONS_tmp = [re.search(r'<a href="(.+?)">(.+?)</a>', ee).group(2) if re.search(r'<a href="(.+?)">(.+?)</a>', ee) is not None else ee for ee in e]
            CONDITIONS_LINK_tmp = ['https://go.drugbank.com' + re.search(r'<a href="(.+?)">(.+?)</a>', ee).group(1) if re.search(r'<a href="(.+?)">(.+?)</a>', ee) is not None else ee for ee in e]
            CONDITIONS.append(';'.join(CONDITIONS_tmp))
            CONDITIONS_LINK.append(';'.join(CONDITIONS_LINK_tmp))
            
        COUNT_ls = [e[4] for e in data['data']]
        
        COUNT = [re.search(r'<a href="(.+?)">(.+?)</a>', e).group(2) for e in COUNT_ls]
        COUNT_LINK = ['https://go.drugbank.com' + re.search(r'<a href="(.+?)">(.+?)</a>', e).group(1) for e in COUNT_ls]

        df = pd.DataFrame({"PAHSE": PHASE,
                           "STATUS": STATUS,
                          "PURPOSE": PURPOSE,
                          "CONDITIONS": CONDITIONS,
                          "CONDITIONS_LINK": CONDITIONS_LINK,
                          "COUNT": COUNT,
                          "COUNT_LINK": COUNT_LINK})
        df_ls.append(df)
        
        start += 50
        if start > data['recordsTotal']:
            break
    
    with open('drugbank_crawler_error_url.txt', 'a') as f:
        f.write('\n'.join(error_url_ls))
        f.write('\n')
    
    df_merge = pd.concat(df_ls)
    return df_merge