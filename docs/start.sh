
nohup inotifywait -m -r . --exclude _build/ --event create,modify | while read changed; do echo $changed; make html>/dev/null; done &

uvicorn --reload --reload-dir _build --port 5000 --host 0.0.0.0 build:app