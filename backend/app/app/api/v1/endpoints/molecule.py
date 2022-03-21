import sys
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


@router.put("/search", response_model=List[schemas.Molecule])
def search_molecules(
    db: Session = Depends(deps.get_db),
    substructure: str = "",
    skip: int = 0,
    limit: int = 100,
    with_dft: bool = True,
    with_xtb: bool = True,
    with_ml: bool = False,
):
    fields = ["dft_data", "xtb_data", "ml_data"]
    bools = [with_dft, with_xtb, with_ml]
    filters = []
    for field, include in zip(fields, bools):
        if include:
            filters.append(f"{field} is not null")
    if filters:
        filters = " or ".join(filters)
    else:
        filters = ""

    limit = limit if limit < 100 else 100
    mol = Chem.MolFromSmiles(substructure)
    if mol is None:
        raise HTTPException(status_code=400, detault="Invalid Smiles Substructure")
    substructure = Chem.MolToSmiles(mol)
    if substructure is None:
        raise HTTPException(status_code=400, default="Invalid Smiles Substructure!")

    sql = text(
        """
        with m as 
            ( select molecule_id,
                     smiles,
                     molecular_weight, 
                     dft_data, 
                     xtb_data,
                     xtb_ni_data, 
                     ml_data 
             from    molecule 
             order by 
                    smiles <-> :substructure 
             offset :offset 
             fetch next :limit rows only ) 
        select m.*,
               array_agg(conformer.conformer_id) as "conformers_id"
        from   m
        left join conformer on (conformer.molecule_id = m.molecule_id)
        group by (m.molecule_id, smiles, molecular_weight,
                  dft_data, xtb_data, xtb_ni_data, ml_data);
        """
    )
    try:
        results = db.execute(
            sql, dict(substructure=substructure, offset=skip, limit=limit)
        ).fetchall()
    except exc.DataError:
        raise HTTPException(status_code=400, detail="Invalid Smiles Substructure!")
    return [schemas.Molecule(**res) for res in results]
