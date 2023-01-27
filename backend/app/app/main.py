from fastapi import FastAPI

import api.v1.api
import api.v2.api

import core.config

app = FastAPI(
    title=core.config.settings.PROJECT_NAME,
    openapi_url="/api/openapi.json",
    redoc=None
)

# versioned API endpoints
app.include_router(api.v1.api.router, prefix="/api/v1")
app.include_router(api.v2.api.router, prefix="/api/v2")

# default API endpoints
app.include_router(api.v2.api.router, prefix="/api")
