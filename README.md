# Simulation Of A Receiving-End Communication System
This project involves five simulations that cover different methods of data recovery for a communication receiver using signal processing techniques. 

The images were formatted to be scanned in columns (with their given unique words for each packet as seen in the bottom of the image) and are modulated using 256-QAM. 
Certain simulations are formateed to have different signal distortions, to simulate real-world data transmission at a receiver. During the receiver simulation, the data 
is received, run through a matched filter and interpolated to handle noise reduction. At this point the data is rotated and placed into a decision matrix to determine if 
the data fits inside a 256-QAM. These can be seen in the plotted "constellations" of the images below, where the blue is the 256-QAM bitmap, and the 'x' marks are the data positions.
If the data does not modulate to the 256-QAM smoothly, the data is placed in feedback systems to estimate phase or timing offset. The plotted error of these estimators is also included in the images 
according to their simulated data's distortion effects.

These simulations were followed with the following criterion (for each we present the results down below):
* Simulation 1: No phase or timing offset through transmission. 
* Simulation 2: Fixed phase offset during transmission, no timing offset
* Simulation 3: Carrier frequency offset, no timing offset
* Simulation 4: No phase offset, fixed timing offset
* Simulation 5: No phase offset, timing clock frequency offset.

## System Results (pigs!)
### Simulation 1
![alt text](https://github.com/17njensen/communication-system-simulator/blob/main/image_results/sim1.PNG)
### Simulation 2
![alt text](https://github.com/17njensen/communication-system-simulator/blob/main/image_results/sim2.PNG)
### Simulation 3
![alt text](https://github.com/17njensen/communication-system-simulator/blob/main/image_results/sim3.PNG)
### Simulation 4
![alt text](https://github.com/17njensen/communication-system-simulator/blob/main/image_results/sim4.PNG)
### Simulation 5
![alt text](https://github.com/17njensen/communication-system-simulator/blob/main/image_results/sim5.PNG)
