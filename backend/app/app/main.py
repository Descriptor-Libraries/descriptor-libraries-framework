
import os

from fastapi import FastAPI

import api.v1.api
import api.v2.api

import core.config

name = os.getenv("API_PREFIX") 

app = FastAPI(
    title=core.config.settings.PROJECT_NAME + " " + name,
    openapi_url=f"/api/{name}/openapi.json",
    docs_url=f"/api/{name}/docs",
    redoc_url=None
)

# Version of the code in this image. The values are set as environment
# variables when the image is built (see backend/Dockerfile); a build that
# does not pass them reports "dev" / "unknown".
@app.get(f"/api/{name}/version", tags=["current"])
def get_version():
    return {
        "version": os.getenv("APP_VERSION", "dev"),
        "git_sha": os.getenv("APP_GIT_SHA", "unknown"),
        "build_date": os.getenv("APP_BUILD_DATE", "unknown"),
    }

# default API endpoints
app.include_router(api.v2.api.router, prefix=f"/api/{name}", tags=["current"])

# versioned API endpoints
app.include_router(api.v2.api.router, prefix=f"/api/{name}/v2", tags=["v2"])
app.include_router(api.v1.api.router, prefix=f"/api/{name}/v1", tags=["v1"])



