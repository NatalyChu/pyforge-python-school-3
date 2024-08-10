from fastapi import FastAPI, status
from fastapi.exceptions import HTTPException
from rdkit import Chem

# from models import User


# [Optional] Upload file with molecules (the choice of format is yours).

app = FastAPI()

molecules = [
    {"mol_id": 1, "name": "CCO"},
    {"mol_id": 2, "name": "c1ccccc1"},
    {"mol_id": 3, "name": "CC(=O)O"},
    {"mol_id": 4, "name": "CC(=O)Oc1ccccc1C(=O)O"},
]


# List all molecules
@app.get("/molecules", tags=["molecules"], summary="List all molecules")
def list_all_molecules(skip: int = 0, limit: int = 10):
    """
    List all molecules
    """
    return molecules[skip : skip + limit]


# Get molecule by identifier
@app.get(
    "/molecules/{mol_id}", tags=["molecules"], summary="Get molecule by identifier"
)
def get_molecule_by_id(mol_id: int):
    """
    Get molecule by identifier
    """
    for mol in molecules:
        if mol["mol_id"] == mol_id:
            return mol
    raise HTTPException(status_code=404, detail="Molecule is not found")


# Add molecule (smiles) and its identifier
@app.post(
    "/add_molecule",
    status_code=status.HTTP_201_CREATED,
    tags=["molecules"],
    response_description="Molecule is added",
)
def add_molecule(molecule: dict):
    for mol in molecules:
        if mol["mol_id"] == molecule["mol_id"]:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Molecule with this ID already exists",
            )
        else:
            molecules.append(molecule)
            return molecule


# Updating a molecule by identifier
@app.put("/molecules/{mol_id}", tags=["molecules"])
def update_molecule(mol_id: int, updated_mol: dict):
    for ind, user in enumerate(molecules):
        if user["mol_id"] == mol_id:
            molecules[ind] = updated_mol
            return updated_mol
    raise HTTPException(status_code=404, detail="Molecule is not found")


# Delete a molecule by identifier
@app.delete(
    "/molecules/{mol_id}",
    tags=["molecules"],
    response_description="Molecule is deleted",
)
def delete_molecule(mol_id: int):
    for ind, mol in enumerate(molecules):
        if mol["mol_id"] == mol_id:
            deleted_mol = molecules.pop(ind)
            return deleted_mol
    raise HTTPException(status_code=404, detail="Molecule is not found")


# Substructure search for all added molecules
@app.post(
    "/molecules/sub_search",
    tags=["molecules"],
    response_description="List of all matching molecules",
)
def find_substructure(sub_str: str):
    """
    Substructure search for all added molecules
    - **sub_str** - a SMILE string to search within all molecules
    """
    match_list = []
    try:
        sub_mol = Chem.MolFromSmiles(sub_str)
        for id, mol in enumerate(molecules):
            smile_mol = Chem.MolFromSmiles(mol["name"])
            if smile_mol.HasSubstructMatch(sub_mol):
                match_list.append(mol["name"])
        return match_list
    except:
        raise HTTPException(
            status_code=422, detail="Provided string cannot be converted into molecule"
        )
