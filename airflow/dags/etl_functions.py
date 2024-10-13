import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from io import BytesIO
import openpyxl
import boto3
from database import get_db
from sqlalchemy.orm import Session

# Function to extract data
def extract_data(**kwargs):
    # Open a session using the existing get_db function
    db: Session = next(get_db())
    # Extract data for the current day
    execution_date = kwargs['ds']  # Get the Airflow execution date
    query = f"SELECT mol_id, name FROM molecules WHERE date = '{execution_date}'"
    result = db.execute(query).fetchall()

    # Convert the result to a DataFrame
    data = pd.DataFrame(result, columns=['mol_id', 'name'])
    
    return data.to_dict()

# Function to transform the data
def transform_data(ti):
    # Fetch extracted data from XComs
    data_dict = ti.xcom_pull(task_ids='extract_data')
    data = pd.DataFrame(data_dict)

    # Calculate molecular properties using RDKit
    smiles = data['name']
    mol_props = []
    
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        mol_weight = Descriptors.MolWt(mol)
        logP = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        h_donors = Descriptors.NumHDonors(mol)
        h_acceptors = Descriptors.NumHAcceptors(mol)
        
        lipinski_pass = (mol_weight < 500 and logP < 5 and h_donors <= 5 and h_acceptors <= 10)
        
        mol_props.append([mol_weight, logP, tpsa, h_donors, h_acceptors, lipinski_pass])
    
    properties_df = pd.DataFrame(mol_props, columns=['Mol_Weight', 'LogP', 'TPSA', 'H_Donors', 'H_Acceptors', 'Lipinski_Pass'])
    
    # Merge with the original data
    transformed_data = pd.concat([data, properties_df], axis=1)
    return transformed_data.to_dict()

# Function to load data to MinIO or S3
def load_data(ti):
    # Fetch transformed data from XComs
    transformed_data_dict = ti.xcom_pull(task_ids='transform_data')
    transformed_data = pd.DataFrame(transformed_data_dict)

    # Save to Excel in memory
    output = BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        transformed_data.to_excel(writer, index=False)

    # Configure  Minio client
    s3_client = boto3.client(
        's3',
        endpoint_url='http://127.0.0.1:9000',
        aws_access_key_id='minioadmin',
        aws_secret_access_key='minioadmin123'
    )
    
    # Upload to MinIO
    s3_client.put_object(
        Bucket='etl-smiles-results',
        Key=f"transformed_data_{ti.execution_date}.xlsx",
        Body=output.getvalue()
    )
