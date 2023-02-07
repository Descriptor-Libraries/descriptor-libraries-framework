from typing import Dict, List, Optional

from pydantic import BaseModel


class Molecule(BaseModel):
    molecule_id: int
    smiles: str
    molecular_weight: Optional[float]
    conformers_id: List[Optional[int]]
    dft_data: Optional[Dict]
    xtb_data: Optional[Dict]
    xtb_ni_data: Optional[Dict]
    ml_data: Optional[Dict]

class MoleculeSimple(BaseModel):
    # Making molecule_id optional for now
    molecule_id: Optional[int]
    umap1: Optional[float]
    umap2: Optional[float]
    pca1: Optional[float]
    pca2: Optional[float]
    pca3: Optional[float]
    pca4: Optional[float]
    pat: Optional[str]
    dist: Optional[float]
    smiles: str
