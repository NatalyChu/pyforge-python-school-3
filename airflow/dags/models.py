from sqlalchemy import Column, Integer, String, Date
from app.database import Base

class Molecule(Base):
    __tablename__ = "molecules"

    mol_id = Column(Integer, primary_key=True, index=True)
    name = Column(String, index=True)
    date = Column(Date, index=True)

    class Config:
        from_attributes = True
