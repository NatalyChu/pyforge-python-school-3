import pytest
from rdkit import Chem
from src.sub_search import substructure_search

# A set of molecules to be used in tests
@pytest.fixture
def molecule_data():
    return [
        "CCO", 
        "c1ccccc1", 
        "CC(=O)O", 
        "CC(=O)Oc1ccccc1C(=O)O"
    ]

# Run the function with different inputs (with given expected outputs) to check how it performs
@pytest.mark.parametrize("substructure, expected", [
    ("c1ccccc1", ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]),
    ("CCO", ["CCO", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]),
    ("CC(=O)O", ["CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]),
    ("N", []),
])
def test_substructure_search(molecule_data, substructure, expected):
    result = substructure_search(molecule_data, substructure)
    expected_mols = [Chem.MolFromSmiles(smiles) for smiles in expected]
    
    assert len(result) == len(expected_mols)
    for res_mol in result:
        assert any(res_mol.HasSubstructMatch(exp_mol) for exp_mol in expected_mols)

# Check that the function raises an exception for an empty or invalid input
@pytest.mark.parametrize("invalid_substructure", [
    None,
    "",
    "invalid_smiles"
])
def test_invalid_substructure_search(molecule_data, invalid_substructure):
    with pytest.raises(Exception):
        substructure_search(molecule_data, invalid_substructure)