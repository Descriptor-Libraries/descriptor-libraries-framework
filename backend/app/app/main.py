from fastapi import FastAPI

import api.v1.api
import core.config

app = FastAPI(
    title=core.config.settings.PROJECT_NAME,
    openapi_url=f'{core.config.settings.API_V1_STR}/openapi.json',
    redoc=None
)

app.include_router(api.v1.api.router, prefix=core.config.settings.API_V1_STR)
