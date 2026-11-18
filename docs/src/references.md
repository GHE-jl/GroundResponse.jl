# References

The models implemented in this package are drawn from the following sources.

## Ground thermal response models

- **Ingersol, L. R.** (1948). Theory of the ground pipe heat source for the heat pump.
  *ASHVE Journal Section, Heating, Piping and Air Conditioning.*
  The infinite line source (ILS) — the exponential-integral solution used by [`ils`](@ref).

- **Carslaw, H. S., & Jaeger, J. C.** (1959). *Conduction of Heat in Solids* (2nd ed.).
  Oxford: Clarendon Press.
  The infinite cylindrical source (ICS) integral evaluated by [`ics`](@ref).

- **Claesson, J., & Javed, S.** (2011). An analytical method to calculate borehole fluid
  temperatures for time-scales from minutes to decades. *ASHRAE Transactions*, 117(PART 2),
  279–288.
  The finite line source (FLS) form implemented in [`fls`](@ref).

- **Pasquier, P., & Lamarche, L.** (2022). Analytic expressions for the moving infinite line source
  model. *Geothermics*, 103, 102413. <https://doi.org/10.1016/j.geothermics.2022.102413>
  The moving infinite line source (MILS), including the early- and late-time series used by
  [`mils`](@ref).

- **Guo, Y., Hu, X., Banks, J., & Liu, W. V.** (2020). Considering buried depth in the moving finite
  line source model for vertical borehole heat exchangers — A new solution. *Energy and Buildings*,
  214, 109859. <https://doi.org/10.1016/j.enbuild.2020.109859>
  The moving finite line source (MFLS) implemented in [`mfls`](@ref).

## Borefield response arrays

- **Cimmino, M., & Bernier, M.** (2014). A semi-analytical method to generate g-functions for
  geothermal bore fields. *International Journal of Heat and Mass Transfer*, 70, 641–650.
  <https://doi.org/10.1016/j.ijheatmasstransfer.2013.11.037>
  The segment-to-segment finite line source of [`fls`](@ref) and the segment array of
  [`gmatrix`](@ref).

- **Cimmino, M.** (2015). The effects of borehole thermal resistances and fluid flow rate on the
  g-functions of geothermal bore fields. *International Journal of Heat and Mass Transfer*, 91,
  1119–1127. <https://doi.org/10.1016/j.ijheatmasstransfer.2015.08.041>
  The BC-IV superposition that uses the segment array of [`gmatrix`](@ref).

- **Dusseault, B., Pasquier, P., & Marcotte, D.** (2018). A block matrix formulation for efficient
  g-function construction. *Renewable Energy*, 121, 249–260.
  <https://doi.org/10.1016/j.renene.2017.12.092>
  The block-matrix formulation whose response array ``G`` is returned by [`gmatrix`](@ref).
