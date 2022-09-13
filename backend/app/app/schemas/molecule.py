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
    smiles: str