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

- **Pasquier, P., Zarrella, A., & Labib, R.** (2018). Application of artificial neural networks to
  near-instant construction of short-term g-functions. *Applied Thermal Engineering*, 143, 910–921.
  <https://doi.org/10.1016/j.applthermaleng.2018.04.078>
  The short-term ANN model (`gST_ANN`), currently not fully integrated.

## Borefield spatial superposition

- **Dusseault, B., Pasquier, P., & Marcotte, D.** (2018). A block matrix formulation for efficient
  g-function construction. *Renewable Energy*, 121, 249–260.
  <https://doi.org/10.1016/j.renene.2017.12.092>
  The block-matrix method behind [`bloc_matrix`](@ref).

- **Nguyen, A., & Pasquier, P.** (2021). A successive flux estimation method for rapid g-function
  construction of small to large-scale ground heat exchanger. *Renewable Energy*, 165, 359–368.
  <https://doi.org/10.1016/j.renene.2020.10.074>
  The successive-flux method behind [`successive_flux`](@ref).
