## coupled-dynamical-systems (CDS)
Implementation of the **Execution** of Coupled Dynamical Systems for ROS &amp; Robot Toolkit. In this specific implementation the **master dynamics** is the position in [x,y,z] and the **slave dynamics** is the orientations [roll,pitch,yaw] of the end-effector of the robot. We use this to learn reaching motions, where a coupling in the position and orientation is present. 

[![Build Status](https://magnum.travis-ci.com/epfl-lasa/coupled-dynamical-systems.svg?token=BqUQb763tsVV4QyzLgBy&branch=master)](https://magnum.travis-ci.com/epfl-lasa/coupled-dynamical-systems)

### Reference:

Shukla, A. and Billard, A. (2011) *Coupled Dynamical System Based Hand-Arm Grasp Planning under Real-Time Perturbations*. In Proceedings of Robotics: Science and Systems VII, Los Angeles CA. Volume 7, Pages 313--320.

### Dependencies:
This package depends on Mathlib which is found in robot-toolkit:
```
$ git clone  https://github.com/epfl-lasa/robot-toolkit.git
```

### Learning CDS:
To learn a CDS with both master/slave dynamics governed by a SEDS model you can use the code in:
```
https://github.com/alpais/cds_model_learning_generic
```

To learn independent SEDS models:
```
https://bitbucket.org/khansari/seds
```
You can find it in SEDS/Extensions/CDSv1, where the GMM parameters of each dynamical systems are stored. These are then fed to CDSExecution class.


### Usage

```C++	
	CDSExecution *cdsRun = new CDSExecution;
	cdsRun->init(masterDyn, slaveDyn, coupling);
	cdsRun->setObjectFrame(object_frame);
	cdsRun->setAttractorFrame(attractor_frame);
	cdsRun->setCurrentEEPose(curr_ee_pose);
	cdsRun->setDT(model_dt);
	cdsRun->setMotionParameters(0.5,1,1,reachingThreshold, masterType, slaveType);
	cdsRun->postInit();
	
	start loop  
		// set current ee pose  
		cdsRun->setCurrentEEPose(curr_ee_pose);  
		// Update DS and get next ee pose   
		cdsRun->getNextEEPose(des_ee_pose); 
	end loop  
```
