from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors


def calculate_descriptors(smiles):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError("Invalid SMILES")

    return {
        "molecular_weight": Descriptors.MolWt(mol),
        "logp": Crippen.MolLogP(mol),
        "molar_refractivity": Crippen.MolMR(mol),
        "heavy_atom_count": mol.GetNumHeavyAtoms(),
        "tpsa": rdMolDescriptors.CalcTPSA(mol)
    }
