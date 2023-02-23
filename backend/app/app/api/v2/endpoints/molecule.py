from multiprocessing.sharedctypes import Value
from re import sub
from typing import List, Optional, Any

from app import schemas
from app.api import deps
from app.db.session import models
from fastapi import APIRouter, Depends, HTTPException
from rdkit import Chem
from sqlalchemy import exc, text
from sqlalchemy.orm import Session

router = APIRouter()

# @router.get("/molecules", response_model=List[schemas.Molecule])
# def get_molecules(db: Session = Depends(deps.get_db)):

#    return "Working"
def valid_smiles(smiles):
    """Check to see if a smile string is valid to represent a molecule.

    Converts the smile string to an rdkit molecule to see if it is valid, then turns it back into a smile string to return an rdkit standardized
    smiles.

    Parameters
    ----------
    smiles : str
             Smile string.
    
    Returns
    -------
    smiles : str
             A smile string generated from an rdkit  molecule.

    Raises
    ------
    HTTPException 
        When an rdkit molecule cannot be created from the smile string.
    HTTPException
        When a smile string cannot be created from an rdkit molecule.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(status_code=400, detail="Invalid Smiles")
    smiles = Chem.MolToSmiles(mol)
    if smiles is None:
        raise HTTPException(status_code=400, detail="Invalid Smiles!")
    
    return smiles

@router.get("/umap", response_model=List[schemas.MoleculeSimple])
def get_molecule_umap(
    limit: int = 1000,
    category: Optional[str] = None,
    show_ml: bool = False,
    db: Session = Depends(deps.get_db),
):

    query_parameters: dict[str, Any] = {"limit": limit}

    query = """
        SELECT molecule_id, smiles, umap[0] AS umap1, umap[1] AS umap2, pat FROM molecule
        """

    # Need to build the query based on the inputs
    # the cleanest way to do this is to consider the
    # possible cases separately.

    # Case 1: no category and no ml (default)
    if not category and not show_ml:
        query += """WHERE dft_data IS NOT NULL OR xtb_data IS NOT NULL
        """

    # Case 2: category and no ml
    if category and not show_ml:
        query += """WHERE pat = :category AND (dft_data IS NOT NULL OR xtb_data IS NOT NULL)
        """
        query_parameters["category"] = category

    # Case 3: category and ml
    if category and show_ml:
        query += """WHERE pat = :category
        """
        query_parameters["category"] = category

    # Case 4: no category and ml
    # no additional arugments needed

    query += """ORDER BY molecule_id 
        FETCH FIRST :limit ROWS ONLY;"""

    sql = text(query)
    stmt = sql.bindparams(**query_parameters)

    results = db.execute(stmt).fetchall()

    return results


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
    # Check if the substructure provided is valid
    substructure = valid_smiles(substructure)

    # The original query is not doing a substructure search at all. It is doing string
    # comparisons between the substructure and the molecule smiles string.
    # added WHERE mol@>:substructure. Leaving order by results in a
    # very slow query, as mentioned in this blog
    # https://depth-first.com/articles/2021/08/11/the-rdkit-postgres-ordered-substructure-search-problem/
    # Adopting their solution works  (set enable_sort=off)
    # However, this isn't the (only) problem.
    # returning the data is very slow.

    # Create mol from smiles in db - select mol_from_smiles('smiles')
    # currently does a substructure search then orders by molecular fingerprint.
    # Timing - without setting enable_sort = off about 34 seconds for 100 molecules
    # with enable_sort off, 0.039 seconds

    sql = text(
        """
        SET LOCAL enable_sort=off;
        select molecule_id, smiles, molecular_weight, umap[0] as umap1, umap[1] as umap2 from molecule 
        where mol@>:substructure
        order by morganbv <%> morganbv_fp(mol_from_smiles(:substructure)) 
        offset :offset 
        limit :limit
        """
    )

    try:
        results = db.execute(
            sql, dict(substructure=substructure, offset=skip, limit=limit)
        ).fetchall()
    except exc.DataError:
        raise HTTPException(status_code=400, detail="Invalid Smiles Substructure!")

    return results

@router.get("/smiles/umap_neighbors/", response_model=List[schemas.MoleculeNeighbors])
def search_umap_neighbors(
    smiles: str,
    skip: int = 1,
    limit: int = 100,
    db: Session = Depends(deps.get_db),
):
    # Check if the smiles provided is valid
    smiles = valid_smiles(smiles)

    sql = text(
        """
        SELECT smiles, ARRAY[umap[0], umap[1]] AS components, molecule_id, 'umap' as type,
        (umap <-> (SELECT umap FROM umap WHERE smiles=:smiles)) as dist
        FROM umap
        ORDER BY dist
        OFFSET :offset 
        LIMIT :limit
        """
    )

    try:
        results = db.execute(
            sql, dict(smiles=smiles, offset=skip, limit=limit)
        ).fetchall()
    except exc.DataError:
        raise HTTPException(status_code=400, detail="No molecule with the smile string provided was found!")

    return results

@router.get("smiles/pca_neighbors/", response_model=List[schemas.MoleculeNeighbors])
def search_pca_neighbors(
    smiles: str,
    components: Optional[str]="1,2,3,4",
    skip: int = 1,
    limit: int = 100,
    db: Session = Depends(deps.get_db),
):
    # Check if the smiles provided is valid
    smiles = valid_smiles(smiles)
    
    query = """"""

    # Need check to see if the components are all integers
    pca_components_list = []
    for i in components.split(","):
        try:
            pca_components_list.append(int(i))
        except ValueError:
            raise HTTPException(status_code=400, detail="Components must be a string of integers separated by commas")

    # Check to see if the pca components requested are valid
    # TODO: This will need to be generalized, since other schemas may have more or less pcas not a maximum of 4...
    pca_components_list.sort()
    if pca_components_list[-1] > 4 or len(pca_components_list) > 4:
        raise HTTPException(status_code=400, detail="Invalid PCA components, there are only 4 available")
        
    query += """SELECT smiles, molecule_id, 'pca' as type, ARRAY["""

    # Look into creating a distance function on SQL to just pass in the parameters

    # Adding on all the PCA columns requested
    for i in pca_components_list:
        if i != pca_components_list[-1]:
            query += f"""pca.pca[{i}],"""
        else:
            query += f"""pca.pca[{i}]] as components,"""

    # Adding on the distance calculation using only PCA columns requested
    for i in pca_components_list:
        if i != pca_components_list[-1]:
            query += f"""+POWER((pca.pca[{i}] - p1.pca[{i}]),2)+"""
        else:
            query += f"""SQRT(POWER((pca.pca[{i}] - p1.pca[{i}]),2)"""

    query += """) as dist
        FROM pca, (SELECT pca FROM pca WHERE smiles=:smiles) as p1
        ORDER BY dist
        OFFSET :offset 
        LIMIT :limit
        """

    sql = text(query)

    try:
        results = db.execute(
            sql, dict(smiles=smiles, offset=skip, limit=limit)
        ).fetchall()
    except exc.DataError:
        raise HTTPException(status_code=400, detail="No molecule with the smile string provided was found!")

    return results

@router.get("/{smiles}/neighbors/{type}/", response_model=List[schemas.MoleculeNeighbors])
def search_neighbors(
    type: str,
    smiles: str,
    components: Optional[str]="1,2,3,4",
    skip: int = 1,
    limit: int = 100,
    db: Session = Depends(deps.get_db),
):
    # Check if the smiles provided is valid
    smiles = valid_smiles(smiles)
    
    query = """"""

    if type == "umap":
        query = """
        SELECT smiles, ARRAY[umap[0], umap[1]] AS components, molecule_id, 'umap' as type,
        (umap <-> (SELECT umap FROM umap WHERE smiles=:smiles)) as dist
        FROM umap
        ORDER BY dist
        OFFSET :offset 
        LIMIT :limit
        """

    if type == "pca":
        # Need check to see if the components are all integers
        pca_components_list = []
        for i in components.split(","):
            try:
                pca_components_list.append(int(i))
            except ValueError:
                raise HTTPException(status_code=400, detail="Components must be a string of integers separated by commas")

        # Check to see if the pca components requested are valid
        # TODO: This will need to be generalized, since other schemas may have more or less pcas not a maximum of 4...
        pca_components_list.sort()
        if pca_components_list[-1] > 4 or len(pca_components_list) > 4:
            raise HTTPException(status_code=400, detail="Invalid PCA components, there are only 4 available")
            
        query += """SELECT smiles, molecule_id, 'pca' as type, ARRAY["""

        # Look into creating a distance function on SQL to just pass in the parameters

        # Adding on all the PCA columns requested
        for i in pca_components_list:
            if i != pca_components_list[-1]:
                query += f"""pca.pca[{i}],"""
            else:
                query += f"""pca.pca[{i}]] as components,"""

        # Adding on the distance calculation using only PCA columns requested
        for i in pca_components_list:
            if i != pca_components_list[-1]:
                query += f"""+POWER((pca.pca[{i}] - p1.pca[{i}]),2)+"""
            else:
                query += f"""SQRT(POWER((pca.pca[{i}] - p1.pca[{i}]),2)"""

        query += """) as dist
            FROM pca, (SELECT pca FROM pca WHERE smiles=:smiles) as p1
            ORDER BY dist
            OFFSET :offset 
            LIMIT :limit
            """

    sql = text(query)

    try:
        results = db.execute(
            sql, dict(smiles=smiles, offset=skip, limit=limit)
        ).fetchall()
    except exc.DataError:
        raise HTTPException(status_code=400, detail="No molecule with the smile string provided was found!")

    return results