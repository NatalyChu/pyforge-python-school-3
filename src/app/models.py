from sqlalchemy import Column, Integer, String
from app.database import Base

class Molecule(Base):
    __tablename__ = "molecules"

    mol_id = Column(Integer, primary_key=True, index=True)
    name = Column(String, index=True)

    class Config:
        from_attributes = True
