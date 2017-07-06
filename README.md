Erdorbit V1.1.0 <br>
Date 06.07.2017 <br>
Author: Thibaut Voirand <br>

See www.erdorbit.com <br>
Plotting orbits in an ECEF frame (Earth centered, Earth fixed) <br>
Inspiration: Quadrature - orbits <br>

Added in this version: <br>

  - refactoring: <br>

      - names modified <br>

      - commentaries reduced <br>

      - functions created for more abstraction: <br>
          - computePositions
          - adaptCoordinatesToCanvasFrame
          - translateCoordinates
          - rotateCoordinates
          - resizeDrawingToFitCanvas
          - drawOrbit
          - getInputParameters

      - code structure improved: <br>
          - first part: variables Declaration
          - second part: math, geometry, and canvas related functions
          - third part: user-interaction related functions
          - fourth part: ajax request to display gallery
