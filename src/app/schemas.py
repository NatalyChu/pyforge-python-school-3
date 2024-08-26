from pydantic import BaseModel

class MoleculeBase(BaseModel):
    name: str

class MoleculeCreate(MoleculeBase):
    pass

class Molecule(MoleculeBase):
    mol_id: int

    class Config:
        orm_mode = True
