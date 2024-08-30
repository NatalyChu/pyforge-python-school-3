import logging
from fastapi import FastAPI, status
from fastapi.exceptions import HTTPException
from rdkit import Chem
from fastapi import FastAPI, Depends, HTTPException
from sqlalchemy.orm import Session
from app import models, schemas
from app.database import engine, get_db
from app.sub_search import substructure_search_db
from typing import Iterator

# Create tables in DB
models.Base.metadata.create_all(bind=engine)

app = FastAPI()

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
@app.post(
    "/molecules/sub_search",
    tags=["molecules"],
    response_description="List of all matching molecules",
    response_model=list[schemas.Molecule]
)
def find_substructure(sub_str: str, db: Session = Depends(get_db)):
    """
    Substructure search for all added molecules
    - sub_str - a SMILE string to search within all molecules
    """
    logger.info(f"Starting substructure search with sub_str: {sub_str}")
    try:
       # Get molecules from the DB
       molecules = db.query(models.Molecule).all()
       logger.info(f"Retrieved {len(molecules)} molecules from the database")

       matching_molecules = substructure_search_db(sub_str, molecules)
       logger.info(f"Found {len(matching_molecules)} matching molecules for sub_str: {sub_str}")

       return matching_molecules

    except Exception as e:
        logger.error(f"Error during substructure search: {str(e)}")
        raise HTTPException(status_code=422, detail="Provided string cannot be converted into molecule")