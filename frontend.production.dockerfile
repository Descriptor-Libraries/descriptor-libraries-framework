FROM node:latest as react-build

WORKDIR /app/
ADD ./frontend /app/
EXPOSE 3000
RUN yarn install
RUN npm run build