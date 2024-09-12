import logging
from fastapi import FastAPI, status
from fastapi.exceptions import HTTPException
from fastapi import FastAPI, Depends, HTTPException
from sqlalchemy.orm import Session
from app import models, schemas
from app.database import engine, get_db
from typing import Iterator
import redis
import json
from tasks import substructure_search_task
from celery.result import AsyncResult

# Create tables in DB
models.Base.metadata.create_all(bind=engine)

app = FastAPI()

# Connect to Redis
redis_client = redis.Redis(host='redis', port=6379, db=0)

# Cache expiration time in seconds
CACHE_EXPIRATION_TIME = 3600

def get_cached_result(key: str):
    result = redis_client.get(key)
    if result:
        return json.loads(result)
    return None

def set_cache(key: str, value: dict, expiration: int = 60):
    redis_client.setex(key, expiration, json.dumps(value))

# Set up logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Add molecule (smiles) and its identifier
@app.post("/molecules/",
    status_code=status.HTTP_201_CREATED,
    tags=["molecules"],
    response_description="Molecule is added",
    response_model=schemas.Molecule)
def create_molecule(molecule: schemas.MoleculeCreate, db: Session = Depends(get_db)):
    logger.info(f"Creating a new molecule with name: {molecule.name}")
    db_molecule = models.Molecule(name=molecule.name)
    db.add(db_molecule)
    db.commit()
    db.refresh(db_molecule)

    logger.info(f"New molecule created with ID: {db_molecule.mol_id} and name: {db_molecule.name}")
    return db_molecule

# List all molecules
def molecule_iterator(db: Session, skip: int, limit: int) -> Iterator[schemas.Molecule]:
    """Iterator returns molecules one by one"""
    query = db.query(models.Molecule).offset(skip).limit(limit)
    for molecule in query:
        yield molecule

@app.get("/molecules/", tags=["molecules"], summary="List all molecules",
          response_model=list[schemas.Molecule])
def read_molecules(skip: int = 0, limit: int = 10, db: Session = Depends(get_db)):
    logger.info(f"Fetching molecules with skip={skip} and limit={limit}")
    iterator = molecule_iterator(db, skip, limit)
    molecules = list(iterator)

    logger.info(f"Fetched {len(molecules)} molecules")
    return molecules


# Get molecule by identifier
@app.get("/molecules/{mol_id}", tags=["molecules"], summary="Get molecule by identifier",
         response_model=schemas.Molecule)
def read_molecule(mol_id: int, db: Session = Depends(get_db)):
    """
    Get molecule by identifier
    """
    logger.info(f"Fetching molecule with ID: {mol_id}")
    db_molecule = db.query(models.Molecule).filter(models.Molecule.mol_id == mol_id).first()
    if db_molecule is None:
        logger.warning(f"Molecule with ID: {mol_id} not found")
        raise HTTPException(status_code=404, detail="Molecule not found")
    
    logger.info(f"Molecule with ID: {mol_id} retrieved successfully")
    return db_molecule

# Delete a molecule by identifier
@app.delete("/molecules/{mol_id}",
            tags=["molecules"],
            response_description="Molecule is deleted", 
            response_model=schemas.Molecule)
def delete_molecule(mol_id: int, db: Session = Depends(get_db)):
    logger.info(f"Attempting to delete molecule with ID: {mol_id}")

    db_molecule = db.query(models.Molecule).filter(models.Molecule.mol_id == mol_id).first()
    if db_molecule is None:
        logger.warning(f"Molecule with ID: {mol_id} not found")
        raise HTTPException(status_code=404, detail="Molecule not found")
    db.delete(db_molecule)
    db.commit()

    logger.info(f"Molecule with ID: {mol_id} deleted successfully")
    return db_molecule

# Updating a molecule by identifier
@app.put("/molecules/{mol_id}", tags=["molecules"],
         response_model=schemas.Molecule)
def update_molecule(mol_id: int, molecule: schemas.MoleculeCreate, db: Session = Depends(get_db)):
    logger.info(f"Attempting to update molecule with ID: {mol_id}")

    db_molecule = db.query(models.Molecule).filter(models.Molecule.mol_id == mol_id).first()
    if db_molecule is None:
        logger.warning(f"Molecule with ID: {mol_id} not found")
        raise HTTPException(status_code=404, detail="Molecule not found")

    db_molecule.name = molecule.name
    db.commit()
    db.refresh(db_molecule)

    logger.info(f"Updating molecule with ID: {mol_id} - new name: {molecule.name}")
    return db_molecule

# Substructure search for all added molecules
@app.post("/molecules/sub_search", tags=["molecules"])
async def start_substructure_search(sub_str: str):
    """
    Starts subsearch task in Celery.
    """
    task = substructure_search_task.delay(sub_str)
    logger.info(f"Started substructure search task with task_id: {task.id}")
    return {"task_id": task.id}

@app.get("/molecules/sub_search/status/{task_id}", tags=["molecules"])
async def get_task_status(task_id: str):
    """
    Gets task status by ID
    """
    task_result = AsyncResult(task_id)
    
    if task_result.state == 'PENDING':
        return {"task_id": task_id, "status": "Pending"}
    elif task_result.state == 'SUCCESS':
        return {"task_id": task_id, "status": "Success", "result": task_result.result}
    elif task_result.state == 'FAILURE':
        return {"task_id": task_id, "status": "Failure", "result": str(task_result.result)}
    
    return {"task_id": task_id, "status": task_result.state}