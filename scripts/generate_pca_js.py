import sys
from tqdm import tqdm
from rdkit import Chem

sys.path.append('../backend/app')
from app.db.session import SessionLocal, models


def get_molecule_id(smiles, session):
    # mol = Chem.MolFromSmiles(smiles)
    # canonical_smiles = Chem.MolToSmiles(mol, canon=True)
    out = (
        session.query(models.molecule)
        .filter(models.molecule.smiles == smiles)
        .one()
    )
    return out.molecule_id

def main(pca_file):
    import ipdb; ipdb.set_trace()
    data = open(pca_file, 'r').read().split('\n')[1:-1]
    dst = open('plot_data.js', 'w')

    session = SessionLocal()

    custom_data = []
    x_list = []
    y_list = []
    for line in tqdm(data):
        _, smiles, x, y = line.split(',')
        molecule_id = get_molecule_id(smiles, session)
        x_list.append(float(x))
        y_list.append(float(y))
        custom_data.append(f'{{ molecule_id: {molecule_id}, smiles: \"{smiles}\" }}')

    dst.write(f"x: {x_list},\n")
    dst.write(f"y: {y_list},\n")
    dst.write(f"customdata: [{', '.join(custom_data)}]\n")
    dst.close()


if __name__ == "__main__":
    main("/Users/tga/Downloads/PL_PC_210126.csv")

