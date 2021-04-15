# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 21:47:16 2021

@author: cmk_8
"""

##################################################################################
######      Query the READ cases - cases Endpoint
##################################################################################
import requests
import json


fields_r1 = [
    "submitter_id",
    "samples.sample_id",
    "samples.sample_type",
    "samples.sample_type_id",
    "samples.submitter_id",
    "samples.is_ffpe",
    "samples.days_to_collection",
    ]

fields_r1 = ",".join(fields_r1)

cases_endpt_r = 'https://api.gdc.cancer.gov/cases'

filters_r1 = {
    "op":"and",
    "content": [
        {
        "op": "in",
        "content":{
            "field": "project.project_id",
            "value": ["TCGA-READ"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "samples.sample_type",
            "value": ["solid tissue normal"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.data_category",
            "value": ["transcriptome profiling"]
            }
        },
        {
        "op":"or",
        "content": [
            {
            "op": "in",
            "content":{
                "field": "files.analysis.workflow_type",
                "value": ["HTSeq - Counts"]
                }
            },
            {
            "op": "in",
            "content":{
                "field": "files.analysis.workflow_type",
                "value": ["HTSeq - FPKM"]
                }
            }],
        }
    ]
}

params_r1 = {
    "filters": json.dumps(filters_r1),
    "fields": fields_r1,
    "format": "JSON",
    "size": "100",
    # "pretty": True
    }

# With a GET request, the filters parameter needs to be converted
# from a dictionary to JSON-formatted string

response_READ_cases = requests.get(cases_endpt_r, params = params_r1)

## build a dictionary with the queried data
dict_READ_cases = response_READ_cases.json()
print(json.dumps((dict_READ_cases), indent=2))


##################################################################################
######      Summarize the data from Query and store data in Dataframe
##################################################################################
import pandas as pd

# Reduce dictionary to directly access data
hits_READ_cases = dict_READ_cases["data"]["hits"]


# Getting rows of hits_samples dictionary and store data in dataframe
rows_READ_cases = []
for data_READ in hits_READ_cases:
   data_row_READ = data_READ["samples"]
   n_READ = data_READ["id"]
   s_READ = data_READ["submitter_id"]
   #f_READ = data_READ[""]

   for row_READ in data_row_READ:
      row_READ["id"] = n_READ
      row_READ["submitter_id"] = s_READ
      rows_READ_cases.append(row_READ)

# Convert to data frame
# df_READ = pd.DataFrame(rows_READ)
# print(df_READ)

## Remove rows with ffpe == True
rows_READ_cases = [i for i in rows_READ_cases if not (i["is_ffpe"] == True)]

## Remove columns with days_to_collection
rows_READ_cases = [i for i in rows_READ_cases if (i["days_to_collection"] == None)]

## Remove columns with sample_type == Blood Derived Normal
rows_READ_cases = [i for i in rows_READ_cases if not (i["sample_type"] == "Blood Derived Normal")]

# Convert to data frame
df_READ_cases = pd.DataFrame(rows_READ_cases)
print(df_READ_cases)

##--------------------------------------------------------------------------------------
"""Here the CMS-Subgroups of the CRC from Synapse are imported
They can later be linked over the submitter_id (=Project) to the querried data"""

## Read in CMS-Subgroups from Synapse
import pandas as pd
df_cms = pd.read_table("cms_labels_public_all.txt")
print(df_cms)
df_cms2 = df_cms
df_cms2.rename(columns = {'sample':'submitter_id'}, inplace = True)

## Add CMS-Subgroups to the READ samples and read it into a .csv file
df_cms_READ = df_READ_cases.merge(df_cms2, on="submitter_id")
df_cms_READ.rename(columns = {'id':'case_id'}, inplace = True)
df_cms_READ
print(df_cms_READ)
## Write the dataframe into a .csv file
#df_cms_READ_csv = df_cms_READ.to_csv("READ.csv",index=True)

#-------------------------------------------------------------------------------------
######################################################################################
#--    Query the FILES endpoint of the READ cases - NORMAL TISSUE ONLY
#######################################################################################
#-------------------------------------------------------------------------------------

import requests
import json

fields_r2 = [
    "file_name",
    "file_id",
    "cases.case_id",
    "cases.submitter_id",
    "cases.samples.sample_type",
    "cases.disease_type",
    "cases.project.project_id",
    "cases.samples.sample_id",
    "cases.samples.sample_type_id",
    "cases.samples.submitter_id",
    "cases.samples.days_to_collection",
    ]

fields_r2 = ",".join(fields_r2)

files_endpt_r = "https://api.gdc.cancer.gov/files"

# This set of filters is nested under an 'and' operator.
filters_r2 = {
        "op":"and",
        "content": [
            {
            "op": "in",
            "content":{
                "field": "cases.project.project_id",
                "value": ["TCGA-READ"]
                }
            },
            {
            "op": "in",
            "content":{
                "field": "cases.samples.sample_type",
                "value": ["solid tissue normal"]
                }
            },
            {
            "op": "in",
            "content":{
                "field": "files.data_category",
                "value": ["transcriptome profiling"]
                }
            },
            {
            "op":"or",
            "content": [
                {
                "op": "in",
                "content":{
                    "field": "files.analysis.workflow_type",
                    "value": ["HTSeq - Counts"]
                    }
                },
                        {
                "op": "in",
                "content":{
                    "field": "files.analysis.workflow_type",
                    "value": ["HTSeq - FPKM"]
                    }
                }],
            }
        ]
}

# A POST is used, so the filter parameters can be passed directly as a Dict object.
params_r2 = {
    "filters": json.dumps(filters_r2),
    "fields": fields_r2,
    "format": "json",
    "size": "1000",
    #"sort": "submitter_id"
    }

# The parameters are passed to 'json' rather than 'params' in this case
response_READ_files_N = requests.get(files_endpt_r, params = params_r2)
dict_READ_files_N = response_READ_files_N.json()
print(response_READ_files_N.content.decode("utf-8"))


##  STORE DATA in DATAFRAME
import pandas as pd

# Build a df with hits-samples
hits_READ_files_N = dict_READ_files_N["data"]["hits"]


# Getting rows of hits_samples dictionary and store data in dataframe
rows_READ_files_N = []
for case in hits_READ_files_N:
    data_case = case["cases"]
    f = case["file_id"]
    n = case["file_name"]
    
    for data in data_case:
       data_row = data["samples"]
       data["file_id"] = f
       data["file_name"] = n
       s = data["submitter_id"]
       c = data["case_id"]
    
       for row in data_row:
          row["file_id"] = f
          row["submitter_id"] = s
          row["file_name"] = n
          row["case_id"] = c
          rows_READ_files_N.append(row)

# Convert to data frame
df_N = pd.DataFrame(rows_READ_files_N)
print(df_N)


## Remove columns with days_to_collection & convert to dataframe
rows_READ_files_N = [i for i in rows_READ_files_N if (i["days_to_collection"] == None)]

df_READ_files_N = pd.DataFrame(rows_READ_files_N)
print(df_READ_files_N)

## Merge the dataframes: df_READ_files_N && df_cms_READ
df_READ_merged_N = pd.merge(df_READ_files_N, df_cms_READ, on=["sample_id", "submitter_id","case_id",
                                                              "sample_type_id", "sample_type"])
df_READ_merged_N




#-------------------------------------------------------------------------------------
######################################################################################
#--    Query the FILES endpoint of the READ cases - PRIMARY Tumor ONLY
#######################################################################################
#-------------------------------------------------------------------------------------

import requests
import json

fields_r3 = [
    "file_name",
    "file_id",
    "cases.case_id",
    "cases.submitter_id",
    "cases.samples.sample_type",
    "cases.disease_type",
    "cases.project.project_id",
    "cases.samples.sample_id",
    "cases.samples.sample_type_id",
    "cases.samples.submitter_id",
    "cases.samples.days_to_collection",
    ]

fields_r3 = ",".join(fields_r3)

files_endpt_r = "https://api.gdc.cancer.gov/files"

# This set of filters is nested under an 'and' operator.
filters_r3 = {
        "op":"and",
        "content": [
            {
            "op": "in",
            "content":{
                "field": "cases.project.project_id",
                "value": ["TCGA-READ"]
                }
            },
            {
            "op": "in",
            "content":{
                "field": "cases.samples.sample_type",
                "value": ["Primary tumor"]
                }
            },
            {
            "op": "in",
            "content":{
                "field": "files.data_category",
                "value": ["transcriptome profiling"]
                }
            },
            {
            "op":"or",
            "content": [
                {
                "op": "in",
                "content":{
                    "field": "files.analysis.workflow_type",
                    "value": ["HTSeq - Counts"]
                    }
                },
                        {
                "op": "in",
                "content":{
                    "field": "files.analysis.workflow_type",
                    "value": ["HTSeq - FPKM"]
                    }
                }],
            }
        ]
}

# A POST is used, so the filter parameters can be passed directly as a Dict object.
params_r3 = {
    "filters": json.dumps(filters_r3),
    "fields": fields_r3,
    "format": "json",
    "size": "1000",
    #"sort": "submitter_id"
    }

# The parameters are passed to 'json' rather than 'params' in this case
response_READ_files_T = requests.get(files_endpt_r, params = params_r3)
dict_READ_files_T = response_READ_files_T.json()
print(response_READ_files_T.content.decode("utf-8"))


##   STORE DATA in DATAFRAME
import pandas as pd

# Build a df with hits-samples
hits_READ_files_T = dict_READ_files_T["data"]["hits"]


# Getting rows of hits_samples dictionary and store data in dataframe
rows_READ_files_T = []
for case in hits_READ_files_T:
    data_case = case["cases"]
    f = case["file_id"]
    n = case["file_name"]
    
    for data in data_case:
       data_row = data["samples"]
       data["file_id"] = f
       data["file_name"] = n
       s = data["submitter_id"]
       c = data["case_id"]
    
       for row in data_row:
          row["file_id"] = f
          row["submitter_id"] = s
          row["file_name"] = n
          row["case_id"] = c
          rows_READ_files_T.append(row)

# Convert to data frame
df_T = pd.DataFrame(rows_READ_files_T)
print(df_T)


## Remove columns with days_to_collection & convert to dataframe
rows_READ_f2_T = [i for i in rows_READ_files_T if (i["days_to_collection"] == None)]

df_READ_files_T = pd.DataFrame(rows_READ_f2_T)
print(df_READ_files_T)

# Reduce Dataframe by merging with df_cms_READ groups
df_READ_merged_T = pd.merge(df_READ_files_T, df_cms_READ, on=["sample_id", "submitter_id","case_id",
                                                              "sample_type_id", "sample_type"])
df_READ_merged_T

#===============================================================================================
# Remove Data entries, which do not have same files for normal tissue and primary tumor
#####################################

df_READ_merged_T_NEW = df_READ_merged_T

def remove_primary_not_in_normal(dataframe):
    exclude_list = []
    i = 0
    new_df = dataframe
    for id_T in dataframe["submitter_id"]:
        if id_T in df_READ_merged_N.values:
            i += 1
            continue
        if id_T not in df_READ_merged_N.values:
            print("Submitter ID " + id_T + "at position ", i, " is not in second Dataframe")
            if id_T not in exclude_list:
                exclude_list.append(id_T)
            new_df = new_df.drop(i)
            i += 1
    print("\nFollowing data was removed from dataframe are: ", exclude_list)
    return new_df

df_READ_merged_T_NEW = remove_primary_not_in_normal(df_READ_merged_T_NEW)

## Merge the dataframes: df_READ_files && df_cms_READ

# df_READ_N_T = pd.merge(df_READ_files_N, df_READ_files_T,how="left", on=["submitter_id", "sample_type"])
# df_READ_N_T
frames = [df_READ_merged_N, df_READ_merged_T_NEW]
df_READ_N_T = pd.concat(frames)


# Rename CMS column
df_READ_N_T.rename(columns = {'CMS_final_network_plus_RFclassifier_in_nonconsensus_samples':'CMS_final'}, inplace = True)
df_READ_N_T.reset_index(inplace=True)


# Count, Look at Distribution of CMS-groupds
cms1=0
cms2=0
cms3=0
cms4=0
nolbl=0
for label in df_READ_N_T["CMS_final"]:
    if label == "CMS1":
        cms1 += 1
    elif label == "CMS2":
        cms2 += 1
    elif label == "CMS3":
        cms3 += 1
    elif label == "CMS4":
        cms4 += 1
    else:
        nolbl+=1
cms_lst= [cms1/4, cms2/4, cms3/4, cms4/4, nolbl/4]
df=pd.DataFrame(cms_lst)
df.index = ["CMS1", "CMS2", "CMS3", "CMS4", "NOLBL"]
df.plot.bar()



#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
## Get the UUID-FileIDs from the merged_df

# # This step populates the download list with the file_ids from the previous query
# file_uuid_list_READ = []
# for file_entry in df_READ_N_T["file_id"]:
#     file_uuid_list_READ.append(file_entry)


# ## Downloading Files

# """An open-access GDC file can be downloaded by appending the file UUID 
# to the data endpoint URL."""

# import requests
# import json
# import re

# ## Downloading a Set of Files based on Filters

# data_endpt_r = "https://api.gdc.cancer.gov/data"

# params_r3 = {"ids": file_uuid_list_READ}

# response_READ_data = requests.post(data_endpt_r, data = json.dumps(params_r3), headers = {"Content-Type": "application/json"})

# response_head_cd_READ = response_READ_data.headers["Content-Disposition"]

# file_name = re.findall("filename=(.+)", response_head_cd_READ)[0]

# with open(file_name, "wb") as output_file:
#     output_file.write(response_READ_data.content)

## ---------------------------------------------------------------------------------
##  Export Phenotype Data as .csv from htseq.count-files for import in DSEq2

## Pre-Clean the Dataframe intended to export "df_READ_N_T"

# df_READ_short = df_READ_N_T.drop(["index", "days_to_collection_x", "days_to_collection_y",
#                                     "is_ffpe", "dataset", "CMS_network",
#                                     "CMS_RFclassifier"], axis=1)

# print (df_READ_short[df_READ_short['file_name'].str.contains('.htseq.counts.gz')])

# df_READ_counts = df_READ_short[df_READ_short['file_name'].str.contains('.htseq.counts.gz')]

# # Sort the df by "file_name" to have same order, as the files are downloaded
# # and imported into the DESeq2 pipeline in R
# df_READ_counts = df_READ_counts.sort_values(by="file_name")

# ## Write the dataframe into a .csv file
# READ_pheno_htseqCounts = df_READ_counts.to_csv("READ_pheno.csv",index=True)




































