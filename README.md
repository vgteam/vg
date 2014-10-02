# vg

## variant graph data structures, interchange, and methods

### building

You'll need the protobuf and jansson development libraries installed on your server.

    sudo apt-get install libprotoc-dev libjansson-dev

Now, obtain the repo and its submodules:

    git clone --recursive https://github.com/ekg/vg.git

Then build with `make`, and test with `test`.
