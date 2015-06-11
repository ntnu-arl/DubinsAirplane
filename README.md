# Dubins Airplane path computation

This repo provides Dubins Airplane model path computation tools to be used within the framework of path/trajectory generation of fixed-wing unmanned aerial vehicles. Dubins airplane is an extension of the classical Dubins car model for the 3D case of an airplane. The specific implementation provided here relies on the formulation presented in:

*Mark Owen, Randal W. Beard and Timothy W. McLain, "Implementing Dubins Airplane Paths on Fixed-Wing UAVs"*

and essentially (as described in this paper) corresponds to a modification of the initial model proposed by Lavalle et al. so that it becomes more consistent with the kinematics of a fixed-wing aircraft. Dubins airplane paths are more complicated than Dubins car paths because of the altitude component. Based on the difference between the altitude of the initial and final configurations, Dubins airplane paths can be classified as low, medium, or high altitude gain. While for medium and high altitude gain there are many different Dubins airplane paths, this implementation selects the pat that maximizes the average altitude throughout the maneuver.

Provided Implementations
-------

**1. Python**

Run the example script "DubinsAirplaneMain.py" in which, one of the 16 supported path-cases may be found (parameter *dubins_case*):

    python DubinsAirplaneMain.py

In the same file, the filename to save the path result (*fname*) can be edited while flags related with plotting (*plot_flag*) and execution verbosity (*verbose_flag*) can be set.

*Example result:*
![enter image description here](https://drive.google.com/file/d/0BwNJOtQI58J0RTRLeENFcVVBSlU/view?usp=sharing)

**Contact**
[Kostas Alexis](mailto:konstantinos.alexis@gmail.com)
