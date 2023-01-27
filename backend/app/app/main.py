from fastapi import FastAPI

import api.v1.api
import api.v2.api

import core.config

app = FastAPI(
    title=core.config.settings.PROJECT_NAME,
    openapi_url="/api/openapi.json",
    redoc=None
)

# default API endpoints
app.include_router(api.v2.api.router, prefix="/api", tags=["current"])

# versioned API endpoints
app.include_router(api.v2.api.router, prefix="/api/v2", tags=["v2"])
app.include_router(api.v1.api.router, prefix="/api/v1", tags=["v1"])



