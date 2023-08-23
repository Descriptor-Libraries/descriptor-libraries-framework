from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.automap import automap_base

from app.core.config import settings


class Models:
    pass


Base = automap_base()
engine = create_engine(str(settings.SQLALCHEMY_DATABASE_URI), pool_pre_ping=True)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base.prepare(engine, reflect=True)

models = Models()
for k in Base.classes.keys():
    setattr(models, k, Base.classes[k])
