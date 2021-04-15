# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 15:39:05 2021

@author: cmk_8
"""

##################################################################################
######      Query the COAD cases - cases Endpoint
##################################################################################
import requests
import json


fields_c1 = [
    "submitter_id",
    "samples.sample_id",
    "samples.sample_type",
    "samples.sample_type_id",
    "samples.submitter_id",
    "samples.is_ffpe",
    "samples.days_to_collection",
    ]

fields_c1 = ",".join(fields_c1)

cases_endpt_c = 'https://api.gdc.cancer.gov/cases'

filters_c1 = {
    "op":"and",
    "content": [
        {
        "op": "in",
        "content":{
            "field": "project.project_id",
            "value": ["TCGA-COAD"]
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

params_c1 = {
    "filters": json.dumps(filters_c1),
    "fields": fields_c1,
    "format": "JSON",
    "size": "1000",
    # "pretty": True
    }

# With a GET request, the filters parameter needs to be converted
# from a dictionary to JSON-formatted string

response_COAD_cases = requests.get(cases_endpt_c, params = params_c1)

## build a dictionary with the queried data
dict_COAD_cases = response_COAD_cases.json()
print(json.dumps((dict_COAD_cases), indent=2))


##  Summarize the data from Query and store data in Dataframe
import pandas as pd

# Reduce dictionary to directly access data
hits_COAD_cases = dict_COAD_cases["data"]["hits"]


# Getting rows of hits_samples dictionary and store data in dataframe
rows_COAD_cases = []
for data_COAD in hits_COAD_cases:
   data_row_COAD = data_COAD["samples"]
   n_COAD = data_COAD["id"]
   s_COAD = data_COAD["submitter_id"]

   for row_COAD in data_row_COAD:
      row_COAD["id"] = n_COAD
      row_COAD["submitter_id"] = s_COAD
      rows_COAD_cases.append(row_COAD)

# Convert to data frame
# df_COAD = pd.DataFrame(rows_COAD_cases)
# print(df_COAD)

## Remove columns with ffpe == True
rows_COAD_cases = [i for i in rows_COAD_cases if not (i["is_ffpe"] == True)]

## Remove columns with days_to_collection
rows_COAD_cases = [i for i in rows_COAD_cases if (i["days_to_collection"] == None)]

## Remove columns with sample_type == Blood Derived Normal
rows_COAD_cases = [i for i in rows_COAD_cases if not (i["sample_type"] == "Blood Derived Normal")]

# Convert to data frame
df_COAD_cases = pd.DataFrame(rows_COAD_cases)
# print(df_COAD_cases)


##--------------------------------------------------------------------------------------
"""Here the CMS-Subgroups of the CRC from Synapse are imported
They can later be linked over the submitter_id (=Project) to the querried data"""

## Read in CMS-Subgroups from Synapse
import pandas as pd
df_cms = pd.read_table("cms_labels_public_all.txt")
# print(df_cms)
df_cms2 = df_cms
df_cms2.rename(columns = {'sample':'submitter_id'}, inplace = True)

## Add CMS-Subgroups to the COAD samples and read it into a .csv file
df_cms_COAD = df_COAD_cases.merge(df_cms2, on="submitter_id")
df_cms_COAD.rename(columns = {'id':'case_id'}, inplace = True)
# print(df_cms_COAD)

## Write the dataframe into a .csv file
#df_cms_COAD_csv = df_cms_COAD.to_csv("COAD.csv",index=True)


#-------------------------------------------------------------------------------------
######################################################################################
#--    Query the FILES endpoint of the READ cases
#######################################################################################
#-------------------------------------------------------------------------------------

import requests
import json

fields_c2 = [
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

fields_c2 = ",".join(fields_c2)

files_endpt_c = "https://api.gdc.cancer.gov/files"

# This set of filters is nested under an 'and' operator.
filters_c2 = {
        "op":"and",
        "content": [
            {
            "op": "in",
            "content":{
                "field": "cases.project.project_id",
                "value": ["TCGA-COAD"]
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
params_c2 = {
    "filters": json.dumps(filters_c2),
    "fields": fields_c2,
    "format": "json",
    "size": "1000",
    #"sort": "submitter_id"
    }

# The parameters are passed to 'json' rather than 'params' in this case
response_COAD_files_N = requests.get(files_endpt_c, params = params_c2)
dict_COAD_files_N = response_COAD_files_N.json()
print(response_COAD_files_N.content.decode("utf-8"))


######      STORE DATA in DATAFRAME
import pandas as pd

# Build a df with hits-samples
hits_COAD_files_N = dict_COAD_files_N["data"]["hits"]


# Getting rows of hits_samples dictionary and store data in dataframe
rows_COAD_files_N = []
for case in hits_COAD_files_N:
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
          rows_COAD_files_N.append(row)

# Convert to data frame
# df = pd.DataFrame(rows_COAD_files_N)
# print(df)


## Remove columns with days_to_collection & convert to dataframe
rows_COAD_f2_N = [i for i in rows_COAD_files_N if (i["days_to_collection"] == None)]

df_COAD_files_N = pd.DataFrame(rows_COAD_f2_N)
print(df_COAD_files_N)

## Merge the dataframes: df_COAD_files_N && df_cms_COAD
df_COAD_merged_N = pd.merge(df_COAD_files_N, df_cms_COAD, on=["sample_id", "submitter_id","case_id",
                                                              "sample_type_id", "sample_type"])
df_COAD_merged_N



#-------------------------------------------------------------------------------------
######################################################################################
#--    Query the FILES endpoint of the COAD cases - PRIMARY Tumor ONLY
#######################################################################################
#-------------------------------------------------------------------------------------

import requests
import json

fields_c3 = [
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

fields_c3 = ",".join(fields_c3)

files_endpt_c = "https://api.gdc.cancer.gov/files"

# This set of filters is nested under an 'and' operator.
filters_c3 = {
        "op":"and",
        "content": [
            {
            "op": "in",
            "content":{
                "field": "cases.project.project_id",
                "value": ["TCGA-COAD"]
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
params_c3 = {
    "filters": json.dumps(filters_c3),
    "fields": fields_c3,
    "format": "json",
    "size": "1000",
    #"sort": "submitter_id"
    }

# The parameters are passed to 'json' rather than 'params' in this case
response_COAD_files_T = requests.get(files_endpt_c, params = params_c3)
dict_COAD_files_T = response_COAD_files_T.json()
print(response_COAD_files_T.content.decode("utf-8"))


######      STORE DATA in DATAFRAME
import pandas as pd

# Build a df with hits-samples
hits_COAD_files_T = dict_COAD_files_T["data"]["hits"]


# Getting rows of hits_samples dictionary and store data in dataframe
rows_COAD_files_T = []
for case in hits_COAD_files_T:
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
          rows_COAD_files_T.append(row)

# Convert to data frame
# df_T = pd.DataFrame(rows_COAD_files_T)
# print(df_T)


## Remove columns with days_to_collection & convert to dataframe
rows_COAD_f2_T = [i for i in rows_COAD_files_T if (i["days_to_collection"] == None)]

df_COAD_files_T = pd.DataFrame(rows_COAD_f2_T)
# print(df_COAD_files_T)

# Reduce Dataframe by merging with df_cms_COAD groups
df_COAD_merged_T = pd.merge(df_COAD_files_T, df_cms_COAD, on=["sample_id", "submitter_id","case_id",
                                                              "sample_type_id", "sample_type"])
df_COAD_merged_T

#===============================================================================================
# Remove Data entries, which do not have same files for normal tissue and primary tumor
#####################################

df_COAD_merged_T_NEW = df_COAD_merged_T

def remove_primary_not_in_normal(dataframe):
    exclude_list = []
    i = 0
    new_df = dataframe
    for id_T in dataframe["submitter_id"]:
        if id_T in df_COAD_merged_N.values:
            i += 1
            continue
        if id_T not in df_COAD_merged_N.values:
            print("Submitter ID " + id_T + " at position ", i, " is not in second Dataframe")
            if id_T not in exclude_list:
                exclude_list.append(id_T)
            new_df = new_df.drop(i)
            i += 1
    print("\nFollowing data was removed from dataframe are: ", exclude_list)
    return new_df

df_COAD_merged_T_NEW = remove_primary_not_in_normal(df_COAD_merged_T_NEW)


## Merge the dataframes: df_READ_files && df_cms_READ
# df_READ_N_T = pd.merge(df_READ_files_N, df_READ_files_T,how="left", on=["submitter_id", "sample_type"])
# df_READ_N_T
frames = [df_COAD_merged_N, df_COAD_merged_T_NEW]
df_COAD_N_T = pd.concat(frames)


# Rename CMS column
df_COAD_N_T.rename(columns = {'CMS_final_network_plus_RFclassifier_in_nonconsensus_samples':'CMS_final'}, inplace = True)
df_COAD_N_T.reset_index(inplace=True)

##############################################################################
# Count, Look at Distribution of CMS-groupds
cms1=0
cms2=0
cms3=0
cms4=0
nolbl=0
for label in df_COAD_N_T["CMS_final"]:
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



## Get the UUID-FileIDs from the merged_df
# This step populates the download list with the file_ids from the previous query
file_uuid_list_COAD = []
for file_entry in df_COAD_N_T["file_id"]:
    file_uuid_list_COAD.append(file_entry)




#-------------------------------------------------------------------------------
###############################################################################
## Downloading a Set of Files based on Filters
##############################################################################
#-------------------------------------------------------------------------------

"""An open-access GDC file can be downloaded by appending the file UUID 
to the data endpoint URL."""

import requests
import json
import re

data_endpt_c = "https://api.gdc.cancer.gov/data"

params_c4 = {"ids": file_uuid_list_COAD}

response_COAD_data = requests.post(data_endpt_c, data = json.dumps(params_c4), headers = {"Content-Type": "application/json"})

response_head_cd_COAD = response_COAD_data.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head_cd_COAD)[0]

with open(file_name, "wb") as output_file:
    output_file.write(response_COAD_data.content)


###---------------------------------------------------------------------------
""" Data will be cleaned and Exported as Phenotype Data as .csv from htseq.count-files 
for import in DSEq2 Pipline """

# Pre-Clean the Dataframe intended to export "df_READ_N_T"

df_COAD_short = df_COAD_N_T.drop(["index", "days_to_collection_x", "days_to_collection_y",
                                    "is_ffpe", "dataset", "CMS_network",
                                    "CMS_RFclassifier"], axis=1)

print (df_COAD_short[df_COAD_short['file_name'].str.contains('.htseq.counts.gz')])

df_COAD_counts = df_COAD_short[df_COAD_short['file_name'].str.contains('.htseq.counts.gz')]

# Sort the df by "file_name" to have same order, as the files are downloaded
# and imported into the DESeq2 pipeline in R
df_COAD_counts = df_COAD_counts.sort_values(by="file_name")

## Write the dataframe into a .csv file
COAD_pheno_htseqCounts = df_COAD_counts.to_csv("COAD_pheno.csv",index=True)


