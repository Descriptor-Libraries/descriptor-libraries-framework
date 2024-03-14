FROM node:alpine as react-build

WORKDIR /app/

COPY ./frontend/package.json ./frontend/yarn.lock /app/

RUN yarn install

COPY ./frontend /app/

RUN npm run build && \
    npm cache clean --force