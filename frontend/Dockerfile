# Node 24 LTS, pinned by version and digest so the base image cannot change underneath us.
FROM node:24-alpine@sha256:61cf3202592df80f42ff3a595636e97039b351d1bd2a1c475f897f328e4870ed AS react-build

WORKDIR /app/

COPY ./frontend/package.json ./frontend/yarn.lock /app/

# Install exactly what yarn.lock records; fail if it is out of date with package.json.
RUN yarn install --frozen-lockfile --non-interactive && \
    yarn cache clean

COPY ./frontend /app/

RUN npm run build && \
    npm cache clean --force

# What this image contains. Pass these at build time; they are recorded as
# image labels and in dist/version.json. They come last so that changing them does not invalidate the
# layers above.
ARG VERSION=dev
ARG GIT_SHA=unknown
ARG BUILD_DATE=unknown

# The page footer reads this file to show which frontend version is running.
RUN printf '{"version":"%s","git_sha":"%s","build_date":"%s"}\n' "${VERSION}" "${GIT_SHA}" "${BUILD_DATE}" > /app/dist/version.json

LABEL org.opencontainers.image.version=${VERSION} \
      org.opencontainers.image.revision=${GIT_SHA} \
      org.opencontainers.image.created=${BUILD_DATE} \
      org.opencontainers.image.source="https://github.com/Descriptor-Libraries/descriptor-libraries-framework" \
      org.opencontainers.image.url="https://github.com/Descriptor-Libraries/descriptor-libraries-framework" \
      org.opencontainers.image.title="descriptor-frontend" \
      org.opencontainers.image.licenses="BSD-3-Clause" \
      org.opencontainers.image.description="Descriptor Libraries frontend (React app served by Express)"
