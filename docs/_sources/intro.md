# D-HydroLogic
This manual provides information on the D-HYDOLogic that [HydroLogic](https://www.hydrologic.nl) 
developed as part of the [HYDROLIB](https://github.com/Deltares/HYDROLIB) and ROI project.

D-HYDOLogic is a set of python tools that can be used the Deltares D-Hydro suite.
Currently, the toolbox consists of two main parts.
- The ROI-Toolbox contains tools to work with the Randstad Overstromingsmodel Instrumentarium (ROI) models.
- The Inundation Toolbox contains tools to visualize output of inundation simulations.
<!---
```{tableofcontents}
```
-->
::::{grid}
:gutter: 3

:::{grid-item-card} ROI
- [ROI - Stand Alone Service](SAS_run_model.ipynb)
- [ROI - Build Model](SAS_build_model.ipynb)
:::

:::{grid-item-card} Inundation Toolbox
- [Reading a his-file](hisreader_usage.ipynb)
- [Reading a clm-file](clm_example_tol.ipynb)
- [Reading a fou-file](fou_example_tol.ipynb)
- [Reading a map-file](map_example_tol.ipynb)
:::

:::{grid-item-card} API-documentation
- [ROI Toolbox](data_structures_api.md)
- [Inundation Toolbox](IT-API-doc.md)
:::
::::