
from typing import Dict, List, Union

from app import schemas
from app.api import deps
from app.db.session import models
from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel
from rdkit import Chem
from sqlalchemy import exc, text
from sqlalchemy.orm import Session

router = APIRouter()

#@router.get("/molecules", response_model=List[schemas.Molecule])
#def get_molecules(db: Session = Depends(deps.get_db)):

#    return "Working"


@router.get("/{molecule_id}", response_model=schemas.Molecule)
def get_a_single_molecule(molecule_id: int, db: Session = Depends(deps.get_db)):
    molecule = (
        db.query(models.molecule)
        .filter(models.molecule.molecule_id == molecule_id)
        .one()
    )

    response = schemas.Molecule(
        molecule_id=molecule.molecule_id,
        smiles=molecule.smiles,
        molecular_weight=molecule.molecular_weight,
        conformers_id=[c.conformer_id for c in molecule.conformer_collection],
        dft_data=molecule.dft_data,
        xtb_data=molecule.xtb_data,
        xtb_ni_data=molecule.xtb_ni_data,
        ml_data=molecule.ml_data,
    )
    return response


@router.get("/search/", response_model=List[schemas.MoleculeSimple])
def search_molecules(
    substructure: str = "",
    skip: int = 0,
    limit: int = 100,
    db: Session = Depends(deps.get_db),
):
    mol = Chem.MolFromSmiles(substructure)
    if mol is None:
        raise HTTPException(status_code=400, detail="Invalid Smiles Substructure")
    substructure = Chem.MolToSmiles(mol)
    if substructure is None:
        raise HTTPException(status_code=400, detail="Invalid Smiles Substructure!")

    # The original query is not doing a substructure search at all. It is doing string 
    # comparisons between the substructure and the molecule smiles string.
    # added WHERE mol@>:substructure. Leaving order by results in a
    # very slow query, as mentioned in this blog 
    # https://depth-first.com/articles/2021/08/11/the-rdkit-postgres-ordered-substructure-search-problem/
    # Adopting their solution works  (set enable_sort=off)
    # However, this isn't the (only) problem. 
    # returning the data is very slow.
    sql = text(
        """
        SET LOCAL enable_sort=off;
        with m as 
            ( select molecule_id,
                     smiles,
                     molecular_weight
             from    molecule 
             where mol@>:substructure
             order by smiles <-> :substructure 
             offset :offset 
             fetch next :limit rows only ) 
        select m.*,
               array_agg(conformer.conformer_id) as "conformers_id"
        from   m
        left join conformer on (conformer.molecule_id = m.molecule_id)
        group by (m.molecule_id, smiles, molecular_weight);
        """
    )
    try:
        results = db.execute(
             sql, dict(substructure=substructure, offset=skip, limit=limit)
         ).fetchall()
    except exc.DataError:
        raise HTTPException(status_code=400, detail="Invalid Smiles Substructure!")

    return results
