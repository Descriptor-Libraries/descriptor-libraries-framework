from fastapi import APIRouter, Depends
from fastapi.responses import PlainTextResponse
from typing import List

from sqlalchemy import text
from sqlalchemy.orm import Session

from app.api import deps
from app.db.session import models
from app import schemas, utils


router = APIRouter()


def get_conformer_and_format(conformer_id, format, db):
    if not format in ['xyz', 'sdf', 'pdb', 'mol']:
        return 'Not implemented yet'

    try:
        coords, elements = db.query(models.conformer.coords,
                                    models.conformer.elements)\
            .filter(models.conformer.conformer_id == conformer_id).one()
    except:
        return None

    xyz = f"{len(coords)}\n\n"
    for atom, (x, y, z) in zip(elements, coords):
        xyz += f"{atom}\t{x}\t{y}\t{z}\n"

    if format == 'xyz':
        return xyz
    
    mol = utils.obabel('xyz', format, xyz)
    return mol


@router.get('/export/{format}/{conformer_id}',
            response_class=PlainTextResponse)
def export_conformer(
    conformer_id: int = 0,
    format: str = 'xyz',
    db: Session = Depends(deps.get_db),
):
    return get_conformer_and_format(conformer_id, format, db)


@router.get('/export/{filename}',
            response_class=PlainTextResponse)
def format_for_ngl(
    filename: str,
    db: Session = Depends(deps.get_db),
):
    conformer_id, format = filename.split('.')
    return get_conformer_and_format(int(conformer_id), format, db)


@router.get('/others_id/{conformer_id}')
def get_other_conformers_id(
    conformer_id: int = 0,
    db: Session = Depends(deps.get_db)
):
    sql = text((
        "select conformer_id from conformer "
        "where molecule_id = ("
        "select molecule_id from conformer where conformer_id = :conformer_id" # modified slightly from true v1
        ");"
    ))
    results = db.execute(sql, conformer_id=conformer_id).fetchall()
    return [r['conformer_id'] for r in results]
