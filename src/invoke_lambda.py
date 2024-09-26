import boto3
import json
import os

def get_aws_credentials():
    # Check if running in GitHub Actions (CI environment)
    if os.getenv('GITHUB_ACTIONS'):
        return {
            'aws_access_key_id': os.getenv('AWS_ACCESS_KEY_ID'),
            'aws_secret_access_key': os.getenv('AWS_SECRET_ACCESS_KEY')
        }
    else:
        # Load credentials from a local file if running locally
        with open('aws_config.json', 'r') as config_file:
            return json.load(config_file)

# Get AWS credentials
aws_creds = get_aws_credentials()

# Initialize a session using the credentials
session = boto3.Session(
    aws_access_key_id=aws_creds['aws_access_key_id'],
    aws_secret_access_key=aws_creds['aws_secret_access_key'],
    region_name='us-east-1'
)

# Initialize the Lambda client
lambda_client = session.client('lambda')

# Define the event you want to send
event = {
    "names": ["Alice", "Bob", "Charlie"]
}

# Invoke the Lambda function
response = lambda_client.invoke(
    FunctionName='HelloStudentFunction',
    InvocationType='RequestResponse',
    Payload=json.dumps(event)
)

# Read the response
response_payload = json.loads(response['Payload'].read())

# Print the response
print(json.dumps(response_payload, indent=4))
