from app.database import get_db
from app.sub_search import substructure_search_db
from app import models, schemas
from celery_app import celery_app

@celery_app.task
def substructure_search_task(sub_str: str):
    db = next(get_db())
    
    # Get molecules from the DB
    molecules = db.query(models.Molecule).all()

    # SQLAlchemy objects to Pydantic model
    molecules_data = []
    for mol in molecules:
        molecule_data = schemas.Molecule.model_validate({
            "mol_id": mol.mol_id,
            "name": mol.name
        }).model_dump()
        molecules_data.append(molecule_data)

    # Search
    matching_molecules = substructure_search_db(sub_str, molecules_data)

    return matching_molecules
