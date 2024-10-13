# Homework 16

# Airflow

Write an Airflow DAG, which would read data from the table in you project and write results to S3 bucket (use minio if no AWS infra is set). DAG should contain at least 3 tasks:

Extract data, which would extract SMILES and related columns from a table for the current day
Transform data, which would transform data, adding column Molecular weight, logP, TPSA, H Donors, H Acceptors and  Lipinski pass properties
Save resulting data as .xlsx file and load it to S3
Schedule your DAG to run a daily basis. 

You might get extra points for:

using the appropriate functionality of Airflow 
Using functions we haven't discuss on a lections (for example, ti properties to get an info about the execution date)
testing your DAGs
keeping clear git history