
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors


def calculate_descriptors(smiles):
    """
    Calculate molecular descriptors using RDKit.
    """

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError("Invalid SMILES structure")

    return {
        "molecular_weight": Descriptors.MolWt(mol),
        "logp": Crippen.MolLogP(mol),
        "molar_refractivity": Crippen.MolMR(mol),
        "heavy_atom_count": mol.GetNumHeavyAtoms(),
        "tpsa": rdMolDescriptors.CalcTPSA(mol)
    }


def analyze_drug_rules(smiles):
    from modules.drug_rules import ghose_filter, egan_rule

    properties = calculate_descriptors(smiles)

    ghose = ghose_filter(
        mw=properties["molecular_weight"],
        logp=properties["logp"],
        atoms=properties["heavy_atom_count"],
        molar_refractivity=properties["molar_refractivity"]
    )

    egan = egan_rule(
        logp=properties["logp"],
        tpsa=properties["tpsa"]
    )

    return {
        "properties": properties,
        "ghose": ghose,
        "egan": egan
    }
