from rdkit import Chem
from rdkit.Chem import Draw


def visualize(molecules: list) -> None:
    # Draw the molecules in a grid
    img = Draw.MolsToGridImage(
        molecules, molsPerRow=2, subImgSize=(200, 200), returnPNG=False
    )
    img.show()

def substructure_search_db(sub_str: str, molecules):
        # Make SMILES str into molecules in RDKit
        substructure = Chem.MolFromSmiles(sub_str)
        if substructure is None:
            raise ValueError("Invalid substructure SMILES")
        
        # Sub structure search
        matching_molecules = []
        for molecule in molecules:
            mol = Chem.MolFromSmiles(molecule.name)
            if mol and mol.HasSubstructMatch(substructure):
                matching_molecules.append(molecule)

        return matching_molecules

def substructure_search(mols: list, mol: str) -> list:
    """takes two arguments: 1) List of molecules as SMILES strings 2) Substructure as SMILES string
    Function returns a list of all molecules from argument 1 that contain substructure from argument 2
    """

    # Handle invalid inputs
    if not mol:
        raise ValueError("Substructure SMILES string is empty or None.")

    # Create a molecules from strings
    smile_mol = Chem.MolFromSmiles(mol)
    molecules = [Chem.MolFromSmiles(smiles) for smiles in mols]

    match_list = []
    for element in molecules:
        # Perform the substructure search
        if element.HasSubstructMatch(smile_mol):
            match_list.append(element)
    
    # visualize(match_list)

    return match_list


assert substructure_search(
    ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1"
)
assert substructure_search(
    ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "CCO"
)
