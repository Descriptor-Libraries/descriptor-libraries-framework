#!/bin/bash

find ./build -type f -exec sed -i 's|/base_url|'"$VITE_BASE_URL"'|g' {} +

npm run start