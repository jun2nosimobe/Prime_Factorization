emcc ecm.cpp -O3 \
  -I /path/to/boost \
  --bind \
  -o ecm_module.js \
  -s WASM=1 \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s SINGLE_FILE=1