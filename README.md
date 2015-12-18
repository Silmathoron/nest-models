# NEST-models (nngt_module) #


## Content ##

Nest module including several model implementations
* Integrate and fire neurons with precise spike timing (alpha conductance or current)
* Usual AEIF or AdExp (Adaptive Exponential Integrate and Fire) neuron
* Grid-precise AEIF (with interpolation to find the exact threshold crossing but with on-grid spikes)
* Precise spike timing AEIF (with interpolation and off-grid spikes)


## Installation ##

NB: you need to keep your "install" nest folder to be able to install modules (compiler needs access to the header files of NEST).

Download and go to the right folder:
`git clone https://github.com/Silmathoron/nest-models.git && cd nest-models`

Then make the install:
        sh bootstrap.sh
        mkdir .build && cd .build
        ../configure --with-nest=/your/path/to/nest-install-folder/bin/nest-config
        make && make install
