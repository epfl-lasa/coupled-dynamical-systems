# coupled-dynamical-systems
Implementation of the Coupled Dynamical Systems for ROS &amp; Robot Toolkit. In this specific implementation the **master dynamics** is the position in [x,y,z] and the **slave dynamics** is the orientations [roll,pitch,yaw] of the end-effector of the robot. We use this to learn reaching motions, where a coupling in the position and orientation is present. 

##Dependencies:
This package depends in Mathlib which is found in robot-toolkit:
```
$ git clone  https://github.com/epfl-lasa/robot-toolkit.git
```
