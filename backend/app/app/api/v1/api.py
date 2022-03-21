from fastapi import APIRouter


from .endpoints import molecule, conformer


router = APIRouter()
router.include_router(molecule.router, tags=['Molecule'], prefix='/molecule')
router.include_router(conformer.router, tags=['Conformer'], prefix='/conformer')
