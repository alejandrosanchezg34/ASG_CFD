## ASG_CFD

This CFD code has been developed by Alejandro Sanchez García as his final thesis for the Master's degree in Aerospace Engineering at Universitat Politècnica de Catalunya (UPC).
The final thesis is called *Study of the flow in a nozzle for rocket engines*, and it will be publicly available soon on [UPCommons](https://upcommons.upc.edu/).

This is a MATLAB code designed to be able to implement a computationally affordable, one-dimensional, compressible model for the flow of a rocket engine nozzle. To that end, three
numerical schemes have been implemented: the MacCormack technique, the Godunov method, and the MUSCL-Hancock scheme.

The purpose of developing such a code is to be able to integrate it into algorithms designed to optimize nozzle and supersonic air inlets, among other geometries, described by a certain
area distribution function $A = A(x)$. This may offer a preliminary design tool to derive unconventional shapes intended to enhance the performance of rockets and other flying vehicles.
