from typing import Dict, List, Optional

from pydantic import BaseModel


class Molecule(BaseModel):
    molecule_id: int
    smiles: str
    molecular_weight: Optional[float] = None
    conformers_id: List[Optional[int]]
    dft_data: Optional[Dict] = None
    xtb_data: Optional[Dict] = None
    xtb_ni_data: Optional[Dict] = None
    ml_data: Optional[Dict] = None

class MoleculeSimple(BaseModel):
    molecule_id: int
    umap1: Optional[float] = None
    umap2: Optional[float] = None
    pat: Optional[str] = None
    smiles: str

class MoleculeComponents(BaseModel):
    type: Optional[str] = None
    molecule_id: int
    pat: Optional[str] = None
    smiles: str
    components: Optional[list[float]] = None

class MoleculeNeighbors(MoleculeComponents):
    dist: Optional[float] = None

class MoleculeIdentifiers(BaseModel):
    smiles: str
    InChI: str
    InChIKey: str
