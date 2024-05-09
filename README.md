# Program to run hrldas.exe for downscaling using ERA5 data

1. Input data ERA5
2. Input data (for land use conditions using geo_em.d?.nc) from WRF preprocessing geogrid
3. Extract necessary variables from geo_em.d03.nc file
4. Create setup file (like initial condition file in WRF)
5. Create boundary condition files (like wrfbdy)
6. Put together and modify nameless.hrldas
7. Run hrldas
