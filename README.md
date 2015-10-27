# coupled-dynamical-systems
Implementation of the Coupled Dynamical Systems for ROS &amp; Robot Toolkit. In this specific implementation the **master dynamics** is the position in [x,y,z] and the **slave dynamics** is the orientations [roll,pitch,yaw] of the end-effector of the robot. We use this to learn reaching motions, where a coupling in the position and orientation is present. 

[![Build Status](https://magnum.travis-ci.com/epfl-lasa/coupled-dynamical-systems.svg?token=BqUQb763tsVV4QyzLgBy&branch=master)](https://magnum.travis-ci.com/epfl-lasa/coupled-dynamical-systems)

###Reference:

Shukla, A. and Billard, A. (2011) *Coupled Dynamical System Based Hand-Arm Grasp Planning under Real-Time Perturbations*. In Proceedings of Robotics: Science and Systems VII, Los Angeles CA. Volume 7, Pages 313--320.

###Dependencies:
This package depends in Mathlib which is found in robot-toolkit:
```
$ git clone  https://github.com/epfl-lasa/robot-toolkit.git
```
