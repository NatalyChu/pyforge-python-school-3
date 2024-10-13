from airflow import DAG
from airflow.operators.python import PythonOperator
from airflow.utils.dates import days_ago
from airflow.utils.dates import timedelta
from etl_functions import *

default_args = {
    'owner': 'natalia',
    'depends_on_past': False,
    'email_on_failure': False,
    'email_on_retry': False,
    'retries': 1,
    'retry_delay': timedelta(minutes=5),
}

with DAG(
    'smiles_etl_pipeline',
    default_args=default_args,
    description='ETL DAG for processing SMILES data',
    schedule_interval=timedelta(days=1),
    start_date=days_ago(1),
    tags=['ETL', 'SMILES'],
) as dag:

    extract_task = PythonOperator(
        task_id='extract_data',
        python_callable=extract_data,
        provide_context=True
    )

    transform_task = PythonOperator(
        task_id='transform_data',
        python_callable=transform_data,
        provide_context=True
    )

    load_task = PythonOperator(
        task_id='load_data',
        python_callable=load_data,
        provide_context=True
    )

    extract_task >> transform_task >> load_task

