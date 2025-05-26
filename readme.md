# how to run:
1) compile the solver in build
2) run main.py
3) check out the results


# cell indexing
All indexing shall be j,i, with j being the vertical coordinate and i being the horizontal coordinate.

## off-integer notation for cell walls
Cells:
 - Begin at 0,0 in the lower left corner
 
Cell walls:
 - Begin at -1/2, 0 for the bottom-most border in the lower left cell
 - Begin at 0, -1/2 for the left-most border in the lower left cell.
 