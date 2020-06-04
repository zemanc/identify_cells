# identify_cells

Algorithm that identifies convective cells and calculates height and width of these cells.

## Criterias
*  Minimum vertical velocity in order to be part of a convective cell
*  Maximum velocity of cell must be higher than 99.99% quantile of positive vertical velocities
*  Height of maximum vertical velocity >= 1250 m above surface and <= 12'500 m
*  
**Remark:** These criterias are quite arbitrary and can be andjusted.

## Use for other fields
The algorithm is quite generic and can be adapted for other fields. But this requires some changes in the code as it's not yet cleaned up and still a prototype.