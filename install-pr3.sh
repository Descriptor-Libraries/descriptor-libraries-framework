echo "Snapshotting filesystem"

echo "Updating source"
git pull origin master

echo "Building frontend"
cd frontend && yarn build

echo "Installing frontend"
cp -fr build/* /zroot/iocage/jails/www/root/usr/local/www/nginx-dist/
cd ..

echo "Installing backend"
cd backend && cp -fr app /zroot/iocage/jails/backend/root/
cd ..
iocage restart backend

echo "Updating DB"
POSTGRES_SERVER=192.168.5.3 POSTGRES_USER=postgres POSTGRES_PASSWD='' POSTGRES_DB=phosphines alembic upgrade head
psql -U postgres -h 192.168.5.3 -d phosphines -c 'REINDEX index molecule_search_mol' ;
