# Use Miniconda3 as the base image
FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /app

# Copy the environment.yml file to the working directory
COPY environment.yml .

# Install dependencies using conda
RUN conda env create -f environment.yml

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

# Copy the application code to the working directory
COPY . .

# Install application dependencies using pip
RUN conda run -n myenv pip install -r requirements.txt

# Expose the port the app runs on
EXPOSE 8000

# Run the FastAPI application
CMD ["conda", "run", "-n", "myenv", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
