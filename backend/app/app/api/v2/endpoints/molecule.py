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

def valid_smiles(smiles):
    """Check to see if a smile string is valid to represent a molecule.

    Converts the smile string to an rdkit molecule to see if it is valid, then turns it back into a smile string to return an rdkit standardized
    smiles.

    Parameters
    ----------
    smiles : str
        Smiles string.
    
    Returns
    -------
    smiles : str
        A smiles string generated from an rdkit  molecule.

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
        SELECT molecule_id, smiles, umap->1 as umap1, umap->2 as umap2, pat FROM molecule
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

@router.get("/{molecule_id}/neighbors/", response_model=List[schemas.MoleculeNeighbors])
def search_neighbors(
    molecule_id: int,
    type: str="pca",
    components: Optional[str]=None,
    skip: int = 1,
    limit: int = 100,
    db: Session = Depends(deps.get_db),
):
    
    type = type.lower()
    
    # Check for valid neighbor type.
    if type not in  ["pca", "umap"]:
        raise HTTPException(status_code=400, detail="Invalid neighbor type.")

    # Set defaults for components
    if type == "pca" and components is None:
        components =  "1,2,3,4"

    # Set defaults for components
    if type == "umap" and components is None:
        components = "1,2"

    # Check correct format
    components_list = []
    for i in components.split(","):
        try:
            components_list.append(int(i))
        except ValueError:
            raise HTTPException(status_code=400, detail="Components must be a string of integers separated by commas")

    # Generalized - get max number of components for type by getting one record
    query = text(f"SELECT cube_dim({type}) FROM molecule WHERE {type} IS NOT NULL LIMIT 1")

    max_dims =  db.execute(query).fetchall()[0][0]

    # Check to see if the pca components requested are valid
    components_list.sort()
    if components_list[-1] > max_dims or len(components_list) > max_dims:
        raise HTTPException(status_code=400, detail=f"Invalid components, there are only {max_dims} available")
    
    query = """"""

    # Creates list of strings of indexing the cube using the `->` operator, ex. ["p2->1", "p2->2", ...]
    cube_indexing = ["p2->" + str(i) for i in range(1, len(components_list)+1)]
    # Creates the string array from the indexing strings. ex. "ARRAY[p2->1, p2->2]"
    array_substitute_one = f'ARRAY[{", ".join(i for i in cube_indexing)}]'
    # Creates the string array for indexing the cube using cube_subset. ex "ARRAY[1, 2, 3, 4]"
    array_substitute_two = f'ARRAY[{", ".join(str(i) for i in components_list)}]'
        
    query = f"""
    SELECT 
        smiles, 
        molecule_id, 
        pat, 
        {array_substitute_one} as components, 
        '{type}' as type, 
        cube_subset(p1.{type}, {array_substitute_two}) <-> p2 as dist
    FROM 
        molecule, 
        (SELECT {type} FROM molecule WHERE molecule_id=:molecule_id) as p1, 
        cube_subset(molecule.{type}, {array_substitute_two}) as p2
    ORDER BY 
        dist
    OFFSET 
        :offset 
    LIMIT 
        :limit
    """

    sql = text(query)

    try:
        results = db.execute(
            sql, dict(molecule_id=molecule_id, offset=skip, limit=limit)
        ).fetchall()
    except exc.DataError:
        raise HTTPException(status_code=400, detail="No molecule with the id provided was found!")

    return results