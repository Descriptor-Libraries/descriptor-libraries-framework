from fastapi import FastAPI

from .api.v1.api import router
from .core.config import settings

app = FastAPI(
    title=settings.PROJECT_NAME,
    openapi_url=f'{settings.API_V1_STR}/openapi.json',
    redoc=None
)

app.include_router(router, prefix=settings.API_V1_STR)
