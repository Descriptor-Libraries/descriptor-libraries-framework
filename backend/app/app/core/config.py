from pydantic import BaseSettings, PostgresDsn, validator

from typing import Any, Optional, Dict


class Settings(BaseSettings):
    PROJECT_NAME: str = 'Phosphines'
    API_V1_STR: str = '/api/v1'

    POSTGRES_SERVER: str
    POSTGRES_USER: str
    POSTGRES_PASSWD: str
    POSTGRES_DB: str
    SQLALCHEMY_DATABASE_URI: Optional[PostgresDsn] = None

    @validator("SQLALCHEMY_DATABASE_URI", pre=True)
    def assemble_db_connection(cls, v: Optional[str], values: Dict[str, Any]) -> Any:
        if isinstance(v, str):
            return v
        return PostgresDsn.build(
            scheme="postgresql",
            user=values.get("POSTGRES_USER"),
            password=values.get("POSTGRES_PASSWD"),
            host=values.get("POSTGRES_SERVER"),
            path=f"/{values.get('POSTGRES_DB') or ''}"
        )

    class Config:
        case_sensitive = True


settings = Settings()
