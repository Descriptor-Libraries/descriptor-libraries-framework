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
    molecule_id: int
    umap1: Optional[float]
    umap2: Optional[float]
    pat: Optional[str]
    smiles: str

class MoleculeNeighbors(BaseModel):
    type: Optional[str]
    molecule_id: int
    dist: Optional[float]
    smiles: str
    components: List[float]


