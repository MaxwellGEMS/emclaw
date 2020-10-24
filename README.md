emclaw
======

Solvers for Maxwell's electrodynamic equations with time-space varying coefficients based on [PyClaw]([https://alals](http://www.clawpack.org/pyclaw/index.html))

# Requierements/Dependencies
- Python (>= 3.7), including recent releases of numpy and scipy 
- [Clawpack](http://www.clawpack.org) (>= 5.7)
- [PETSc](https://gitlab.com/petsc/petsc) and [PETSc4py](https://gitlab.com/petsc/petsc4py)
- ParaView, YT, or Avizo to visualize 3D data (experimental)

# quick _setup_

1. clone the repository
2. compile the Riemann sources
   ```
   export EMCLAW_PATH=path_to_emclaw_clone
   cd $EMCLAW_PATH\emclaw\riemann
   python setup.py build_ext -i
   ```
3. Add emclaw to an environment as a development module (e.g. `emclaw_env`)
   ```
    conda develop -n emclaw_env $EMCLAW_PATH
    ```
4. Test
    ```
    conda activate emclaw_env
    python -c "import emclaw; print(emclaw.__path__[0])"
    >> EMCLAW_PATH
    ```

# Documentation
[Read the docs](https://emclaw.readthedocs.io/en/latest/)

# Questions/Contact
damian.sanroman@kaust.edu.sa


