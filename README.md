Simulation of Simply Supported Beam
====================================
# Simulation of Simply Supported Beam
Description
------------
- This package carry out output based Modal Analysis for an array of sensor network on simply supported beam [link](https://github.com/tahsin-nishat/modal/edit/main/README.md)

Expected output
----------------
1. Go to datagen.py file to generate acceleration response time history based on a simlpy supported beam geometry and material properties.
2. User can define whether the input load is stochastic in nature or impact or both. Based on the load input, acceleration time history varies. User can also predefine the noise level. The expected signal and the power spectral for stochastic load should look like below:
   
   ![acceleration](https://github.com/tahsin-nishat/modal/blob/main/results/Acceleration_at_Node_2.png)

   ![psd](https://github.com/tahsin-nishat/modal/blob/main/results/PSD_of_Displacement_at_last_node.png)

  For an impact load input, user need to define load case (i.e., `load=1` or `load=2`). To visualize the impact load, the following is the reference impact at node 3
  
  ![impact](https://github.com/tahsin-nishat/modal/blob/main/results/Impact_load.png)
    
  The expected signal and the power spectral for impact load should look like below:
  
  ![acceleration](https://github.com/tahsin-nishat/modal/blob/main/results/Acceleration_at_Node_3.png)
    
  ![psd](https://github.com/tahsin-nishat/modal/blob/main/results/PSD_of_impactacceleration_at_last_node.png)
  
2. Go to findrefmod.py to conduct output-based covariance driven stochastic system identification. The file is setup in such a way that user needs one reference response generated from datagen.py and an estimated response generated in other way or just skipping some sensor location. Then one can compare two different modeshapes to see how they relate. Following are the references. The first one is a stabilization diagram where user has to manually pick the poles based on observation and then one can see two different mode shape comparison.

      ![stabilization](https://github.com/tahsin-nishat/modal/blob/main/results/Stabilization_Diagram.png)
   
     ![modeshape](https://github.com/tahsin-nishat/modal/blob/main/results/modeshape_com.png)
