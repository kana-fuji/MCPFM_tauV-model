# MCPFM_tauV_model


The `MCPFM_tauV_model` library implements the multicellular phase-field model.

## Overview

### The multicellular phase-field model


## Usage

```bash

# Compile the code
g++ -I. -std=c++11 -O3 -DSFMT_MEXP=19937 -o run SFMT.c mcpf_2d_us.C -lm

# Run the simulation from 4 cells: 
./run.sh test_n4_tauV50_xi016_rs070 4 50 0.16 0.70

```

- The first argument "test_n4_tauV50_xi016_rs070" represents the parameter name.
- The second argument "4" is the initial cell number.
- The third argument "50" is the growth time scale.
- The fourth argument "0.16" is the lumen pressure.
- The fifth argument "0.70" is the initial lumen size.


## Dependencies

- C++


## Authors

* Kana Fuji - The University of Tokyo
* fuji@g.ecc.u-tokyo.ac.jp


## License
The 'MCPFM_tauV_model' is licensed under the [MIT license](https://en.wikipedia.org/wiki/MIT_License).
