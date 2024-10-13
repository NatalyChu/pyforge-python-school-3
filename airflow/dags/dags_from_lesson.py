import pandas as pd
import pendulum
import requests
from airflow import DAG
from airflow.models import Variable
from airflow.operators.empty import EmptyOperator
from airflow.operators.python import PythonOperator
from airflow.providers.postgres.hooks.postgres import PostgresHook


def download(ti):
    biogrid_url = Variable.get('biogrid_url')
    response = requests.get(
        biogrid_url.format(version='4.4.200'),
        params={'downloadformat': 'zip'}
    )

    if response.status_code == 200:
        local_file_name = 'biogrid.tab3.zip'
        with open(local_file_name, 'wb') as f:
            f.write(response.content)
        ti.xcom_push(
            key='local_file_name',
            value=local_file_name
        )
    else:
        raise Exception('No data was fetched')


def transform_and_ingest(ti):
    local_file_name = ti.xcom_pull(
        task_ids='download_and_transform',
        key='local_file_name'
    )
    df = pd.read_csv(local_file_name, delimiter='\t', compression='zip')

    df = df.rename(
        lambda column_name: column_name.lower().replace(' ', '_').replace('#', '_').strip('_'),
        axis='columns'
    )

    df = df


    df['version'] = '4.4.200'

    postgres_hook = PostgresHook(postgres_conn_id='postgres_local')
    engine = postgres_hook.get_sqlalchemy_engine()
    df.to_sql('biogrid_data', engine, if_exists='replace')


with DAG(
    dag_id='biogrid_library_loading_dag',
    start_date=pendulum.today(),
    schedule=None,
    tags=['python_school', 'biogrid']
) as dag:
    start_op = EmptyOperator(
        task_id='start'
    )
    download_and_transform_op = PythonOperator(
        task_id='download_and_transform',
        python_callable=download
    )
    ingest_op = PythonOperator(
        task_id='ingest',
        python_callable=transform_and_ingest
    )
    finish_op = EmptyOperator(
        task_id='finish'
    )

    start_op >> download_and_transform_op >> ingest_op
    ingest_op >> finish_op