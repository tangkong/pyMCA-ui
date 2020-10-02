SSRL pyMCA viewer
=================


Usage
-----
From outside the scripts folder, run a script from the scripts folder

To run the standard viewer:
```
./scripts/run_me.sh
```

Different `.sh` files run different versions of the UI, which could involve different fitting protocols, etc.  

You may have to `chmod +x *.sh` the files in the scripts folder

Data Files
----------
Currently there are two data files you can open from the MCA viewer, I do not know what the difference is.


Demo IOC
--------

This demo IOC uses pcaspy to simulate an EPICS IOC.

Open a terminal and do:

python demo/testing-ioc

In another terminal:

pydm -m DEVICE=XPRESS_DEMO SSRL_MCA.py
