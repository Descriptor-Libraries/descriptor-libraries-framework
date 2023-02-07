from starlette.applications import Starlette
from starlette.routing import Mount
from starlette.staticfiles import StaticFiles


routes = [
    Mount('/docs/', app=StaticFiles(directory='_build/html', html=True), name="documentation"),

]
 
app = Starlette(routes=routes)