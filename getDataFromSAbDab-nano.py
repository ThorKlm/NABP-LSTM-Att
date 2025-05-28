#This script aims to get the data_prev of Nanobodies and their binding antigens from SAbDab-nano database to 3 csv folders
# each folder contacin number of csvs based on fvs number 

import requests
import csv  
import pandas as pd
from bs4 import BeautifulSoup
import os
import re

def getAllRecords():
    print("Calling All structures page ... ")
    url = requests.get(baseUrl+allRecordsEndPoint)
    htmltext = url.text
    data = pd.read_html(htmltext)[0]
    print(str(data.count()["PDB"]) + " PDBs is loaded now, lets save them....")
    return data

def writeDataToCSV(recordType, data):
    os.makedirs(recordType, exist_ok=True)
    fvsNumber = int(data["basic_data"]["Number of Fvs"])
    empty_6_items = ["","","","","",""]
    empty_3_items = ["","",""]
    file_dir = recordType+'/'+str(fvsNumber)+'.csv' 

    if recordType == "nano":
        headers = ['pdb', 'Hchain', 'Hchain_sequence', 'Heavy subgroup',
                'Species', 'In complex?','scFv?','Has constant domain?',
                'CDRH1', 'CDRH2','CDRH3',
                'Antigen chain_1', 'Antigen chain_2', 'Antigen chain_3',
                'antigen_type','antigen_name','Antigen species', 
                'Antigen sequence_1', 'Antigen sequence_2','Antigen sequence_3',
                'HL', 'HC1', 'HC2','LC1','LC2','dc','Method','Resolution','Has constant region'
                ]
    else:       
        headers = ['pdb', 'Hchain', 'Lchain', 'Hchain_sequence','Lchain_sequence', 'Heavy subgroup',
                'Light subgroup', 'Species', 'In complex?','scFv?','Has constant domain?',
                'CDRH1', 'CDRH2','CDRH3','CDRL1', 'CDRL2','CDRL3',
                'Antigen chain_1', 'Antigen chain_2', 'Antigen chain_3',
                'antigen_type','antigen_name','Antigen species', 
                'Antigen sequence_1', 'Antigen sequence_2','Antigen sequence_3',
                'HL', 'HC1', 'HC2','LC1','LC2','dc','Method','Resolution','Light chain type','Has constant region'
        ]

    with open(file_dir, 'a', encoding='UTF8', newline='') as b:
        wr = csv.writer(b)
        is_empty = os.stat(file_dir).st_size == 0
        if is_empty:
            wr.writerow(headers)

    for i in range(fvsNumber):
        cdrSec = data["fv_"+str(i)].get("CDR Sequences (chothia definition)")
        antigenDetailsSec = data["fv_"+str(i)].get("Antigen Details")
        orientationSec = data["fv_"+str(i)].get("Orientation Angles (from ABangle)")
        if recordType == "nano": 
            rowData = [
                data["basic_data"].get("PDB"), 
                data["fv_"+str(i)]["Fv Details"].get("Heavy chain"), 
                data["fv_"+str(i)].get("heavyChainSeq"), 
                data["fv_"+str(i)]["Fv Details"].get("Heavy subgroup"), 
                data["fv_"+str(i)]["Fv Details"].get("Species"), 
                data["fv_"+str(i)]["Fv Details"].get("In complex?"), 
                data["fv_"+str(i)]["Fv Details"].get("scFv?"), 
                data["fv_"+str(i)]["Fv Details"].get("Has constant domain?"), 
            ]
        else:
            rowData = [
                data["basic_data"].get("PDB"), 
                data["fv_"+str(i)]["Fv Details"].get("Heavy chain"), 
                data["fv_"+str(i)]["Fv Details"].get("Light chain"), 
                data["fv_"+str(i)].get("heavyChainSeq"), 
                data["fv_"+str(i)].get("lightChainSeq"), 
                data["fv_"+str(i)]["Fv Details"].get("Heavy subgroup"), 
                data["fv_"+str(i)]["Fv Details"].get("Light subgroup"), 
                data["fv_"+str(i)]["Fv Details"].get("Species"), 
                data["fv_"+str(i)]["Fv Details"].get("In complex?"), 
                data["fv_"+str(i)]["Fv Details"].get("scFv?"), 
                data["fv_"+str(i)]["Fv Details"].get("Has constant domain?"), 
            ]
        antiGenChain1 = antiGenChain2 = antiGenChain3 = ""
        antiGenSeq1 = antiGenSeq2 = antiGenSeq3 = ""
        if antigenDetailsSec != None:
            antiGenString =  data["fv_"+str(i)]["Antigen Details"]["Antigen chains"]
            if antiGenString != None:
                antiGenChains =  antiGenString.split(",")
                if int(len(antiGenChains)) >= 1:
                    antiGenChain1 = antiGenChains[0]
                if int(len(antiGenChains)) >= 2:
                    antiGenChain2 = antiGenChains[1]
                if int(len(antiGenChains)) >= 3:
                    antiGenChain3 = antiGenChains[2]  

            antiGenstring = data["fv_"+str(i)]["Antigen Details"]["Antigen sequence"]
            if antiGenstring != None:
                antiGenSeqs =  antiGenstring.split("/")
                if int(len(antiGenSeqs)) >= 1:
                    antiGenSeq1 = antiGenSeqs[0]
                if int(len(antiGenSeqs)) >= 2:
                    antiGenSeq2 = antiGenSeqs[1]
                if int(len(antiGenSeqs)) >= 3:
                    antiGenSeq3 = antiGenSeqs[2]  

        if cdrSec != None and recordType != "nano":
            rowData.append(data["fv_"+str(i)]["CDR Sequences (chothia definition)"].get("CDRH1")) 
            rowData.append(data["fv_"+str(i)]["CDR Sequences (chothia definition)"].get("CDRH2")) 
            rowData.append(data["fv_"+str(i)]["CDR Sequences (chothia definition)"].get("CDRH3")) 
            rowData.append(data["fv_"+str(i)]["CDR Sequences (chothia definition)"].get("CDRL1")) 
            rowData.append(data["fv_"+str(i)]["CDR Sequences (chothia definition)"].get("CDRL2")) 
            rowData.append(data["fv_"+str(i)]["CDR Sequences (chothia definition)"].get("CDRL3"))
        if cdrSec != None and recordType == "nano":
            rowData.append(data["fv_"+str(i)]["CDR Sequences (chothia definition)"].get("CDRH1")) 
            rowData.append(data["fv_"+str(i)]["CDR Sequences (chothia definition)"].get("CDRH2")) 
            rowData.append(data["fv_"+str(i)]["CDR Sequences (chothia definition)"].get("CDRH3")) 

        if  cdrSec == None and recordType == "nano":
            rowData.extend(empty_3_items)
        if  cdrSec == None and recordType != "nano":
            rowData.extend(empty_6_items) 

        rowData.append(antiGenChain1)
        rowData.append(antiGenChain2)
        rowData.append(antiGenChain3)

        if antigenDetailsSec != None:
            rowData.append(data["fv_"+str(i)]["Antigen Details"].get("Antigen type"))
            rowData.append(data["fv_"+str(i)]["Antigen Details"].get("Antigen name"))
            rowData.append(data["fv_"+str(i)]["Antigen Details"].get("Antigen species"))
        else:
            rowData.extend(empty_3_items)

        rowData.append(antiGenSeq1)
        rowData.append(antiGenSeq2)
        rowData.append(antiGenSeq3)

        if orientationSec != None:
            rowData.append(removeUnit(data["fv_"+str(i)]["Orientation Angles (from ABangle)"].get("HL")))
            rowData.append(removeUnit(data["fv_"+str(i)]["Orientation Angles (from ABangle)"].get("HC1")))
            rowData.append(removeUnit(data["fv_"+str(i)]["Orientation Angles (from ABangle)"].get("HC2")))
            rowData.append(removeUnit(data["fv_"+str(i)]["Orientation Angles (from ABangle)"].get("LC1")))
            rowData.append(removeUnit(data["fv_"+str(i)]["Orientation Angles (from ABangle)"].get("LC2")))
            rowData.append(removeUnit(data["fv_"+str(i)]["Orientation Angles (from ABangle)"].get("dc")))
        else:
            rowData.extend(empty_6_items)  

        rowData.append(data["basic_data"].get("Method")) 

        rowData.append(removeUnit(data["basic_data"].get("Resolution"))) 
        
        if recordType != "nano":
            rowData.append(data["basic_data"].get("Light chain type")) 

        rowData.append(data["basic_data"].get("Has constant region"))
            
        # print(rowData)
        with open(file_dir, 'a', encoding='UTF8', newline='') as body:
            writer = csv.writer(body)
            writer.writerow(rowData)

def getFvData(pdbFvDiv):
    # tables with table-alignment classes are heavy and light chain seqs
    fvData = {}
    seqTableNumber = 0
    for table in pdbFvDiv.find_all('table', class_='table-results'):
        tableData = {}
        for i, row in enumerate(table.find_all('tr')):
            if i == 0:
                header = row.find('th').text.strip()
            else:
                tds = [row.findAll('td')]
                for td in tds :
                    tableData[td[0].string] = td[1].string
            fvData[header] = tableData

    for table in pdbFvDiv.find_all('table', class_='table-alignment'):
        for i, row in enumerate(table.find_all('tr')):
            if i == 0:
                # skip headers row
                continue
            else:
                if seqTableNumber == 0:
                    heavyChainSeq = row.find_all("td")
                    heavyChainSeqData = ""
                    for chain in heavyChainSeq:
                        heavyChainSeqData += chain.text
                    fvData["heavyChainSeq"] = heavyChainSeqData
                else:
                    lightChainSeq = row.find_all("td")
                    lightChainSeqData = ""
                    for chain in lightChainSeq:
                        lightChainSeqData += chain.text
                    fvData["lightChainSeq"] = lightChainSeqData
        seqTableNumber += 1
    return fvData

def getRecordType(data):
    heavyChainCounter = 0
    lightChainCounter = 0
    # if all fvs have heavy chain only and light chain is always none then it's nano
    # if all fvs have both heavy and light chains then anti
    # else combined 
    i = -1
    for fvs in data.items():
        i+=1
        if i == 0:
            #  skip basic_data
            continue
        for fv in fvs:
            if isinstance(fv, dict):
                if fv['Fv Details']['Heavy chain'] != None:
                    heavyChainCounter+=1
                if fv['Fv Details']['Light chain'] != None:
                    lightChainCounter+=1

    if heavyChainCounter == lightChainCounter == int(data["basic_data"]["Number of Fvs"]):
        global antiBodiesCounter
        antiBodiesCounter+=1
        return "anti"
    elif  heavyChainCounter == int(data["basic_data"]["Number of Fvs"]) and lightChainCounter == 0:
        global nanoBodiesCounter
        nanoBodiesCounter+=1
        return "nano"
    else:
        global combinedBodiesCounter
        combinedBodiesCounter+=1
        return "combined"

def removeUnit(value):
    if value == None:
        return value
    else:
        newValue = re.findall(r"[-+]?(?:\d*\.*\d+)", value)
        return newValue[0]

def is_pdb_processed(recordType, fvsNumber, pdb_code):
    file_path = os.path.join(recordType, f"{fvsNumber}.csv")
    if not os.path.exists(file_path):
        return False
    with open(file_path, 'r', encoding='utf8') as f:
        return any(pdb_code in line for line in f)

def get_processed_pdbs():
    processed = set()
    for folder in ['nano', 'anti', 'combined']:
        if not os.path.exists(folder):
            continue
        for filename in os.listdir(folder):
            if filename.endswith('.csv'):
                with open(os.path.join(folder, filename), encoding='utf8') as f:
                    next(f)  # skip header
                    for line in f:
                        pdb_code = line.split(',')[0].strip()
                        processed.add(pdb_code)
    return processed

os.makedirs('./data/nano', exist_ok=True)
os.makedirs('./data/anti', exist_ok=True)
os.makedirs('./data/combined', exist_ok=True)

baseUrl = "https://opig.stats.ox.ac.uk/"
allRecordsEndPoint = "webapps/sabdab-sabpred/sabdab/nano/?all=true"
antiBodiesCounter = 0
nanoBodiesCounter = 0
combinedBodiesCounter = 0
failedBodiesCounter = 0  
recordsCount = 0
allRecords = getAllRecords()
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By

options = webdriver.ChromeOptions()
options.add_argument("--headless")  # Optional: run Chrome in headless mode
processed_pdbs = get_processed_pdbs()

driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=options)
# WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.ID, f"collapse_{fvIndex}")))
os.makedirs("errors_html", exist_ok=True)
from tqdm import tqdm

for index, row in tqdm(allRecords.iterrows(), total=len(allRecords), desc="Processing PDBs"):
    pdb_code = row['PDB']
    if pdb_code in processed_pdbs:
        print(f"{pdb_code} already processed. Skipping.")
        continue
    try:
        recordsCount+=1
        print("Record number "+str(recordsCount) + ": Calling  "+row['PDB']+" page to be saved ... ")
        pdbEndpoint = "webapps/sabdab-sabpred/sabdab/structureviewer/?pdb="+row['PDB']
        # https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/structureviewer/?pdb=8s2t
        # webapps/sabdab-sabpred/sabdab/structureviewer/?pdb="+row['PDB']
        # pdbEndpoint = "webapps/sabdab-sabpred/sabdab/structureviewer/?pdb=1g9e"

        options = webdriver.ChromeOptions()
        options.add_argument("--headless")
        driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=options)
        # WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.ID, f"collapse_{fvIndex}")))

        # driver.get(baseUrl + pdbEndpoint)
        full_url = baseUrl + pdbEndpoint
        print(f"Loading URL: {full_url}")
        driver.get(full_url)
        html = driver.page_source
        soup = BeautifulSoup(html, "lxml")
        # get basic details table 
        pdbDetailsTableDiv = soup.find(id="details")
        pdbDetailsTableRows = pdbDetailsTableDiv.find('table').find_all('tr')
        print("pdbDetailsTableRows", pdbDetailsTableRows)
        wholeData = {}
        tds = [row.findAll('td') for row in pdbDetailsTableRows]
        wholeData["basic_data"] = { td[0].string: td[1].string for td in tds }
        # loop on fvs to get its data_prev
        for fvIndex in range(int(wholeData["basic_data"]["Number of Fvs"])):
            #collapse_0 is the first fv, collapse_1 is the second fv and so on
            WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.ID, f"collapse_{fvIndex}")))
            # collapseDivId = "collapse_"+str(fvIndex)
            collapseDivId = f"collapse_{fvIndex}"
            pdbFvDiv = soup.find(id=collapseDivId)
            if pdbFvDiv is None:
                raise ValueError(f"Missing Fv div for ID: {collapseDivId}")
            # pdbFvDiv = soup.find(id=collapseDivId)

            fvData = getFvData(pdbFvDiv)
            wholeData["fv_"+str(fvIndex)] = fvData

        recordType = getRecordType(wholeData)
        fvsNumber = int(wholeData["basic_data"]["Number of Fvs"])
        pdb_code = wholeData["basic_data"]["PDB"]

        if is_pdb_processed(recordType, fvsNumber, pdb_code):
            print(f"{pdb_code} already processed. Skipping.")
            continue

        writeDataToCSV(recordType, wholeData)
    except Exception as e:
        with open(f"errors_html/{row['PDB']}.html", "w", encoding="utf-8") as f:
            f.write(driver.page_source)
        with open('errors.txt', 'a') as log:
            log.write(row['PDB']+ " falied to be loaded > "+str(e))
            failedBodiesCounter+=1
            log.write('\n')
            print(row['PDB']+" Failed, Fail number: "+str(failedBodiesCounter)+ " till now")
            continue
print("I am done with follwing scraps")
print("Nanobodies counter is " + str(nanoBodiesCounter))
print("Antibodies counter is " + str(antiBodiesCounter))
print("Combined counter is " + str(combinedBodiesCounter))
print("Total fails counter is " + str(failedBodiesCounter))

