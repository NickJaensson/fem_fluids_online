# fem_fluids_online
This repo contains the source code to build a Jupyter book for
*FEM for fluids*.

The book can be accessed here: 
https://nickjaensson.github.io/fem_fluids_online


## Building the book
Install required packages, by running
```bash
python -m pip install -r requirements.txt
```
Compile the book by running in the root directory of the project
```bash
jupyter-book build .
```
The book will be compiled in html format in the `_build` directory.
