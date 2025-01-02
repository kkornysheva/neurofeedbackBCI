# **Neurofeedback**

## Folder Structure 
```
│projectdir          <- Project's main folder
│
├── examples         <- Example scripts
│   ├── brainflow    <- Brainflow is a library to obtain and analyse EEG data, often 
│   │                   recommended for OpenBCI, but I didn't use it yet
│   ├── pylsl        <- Python interface to the LabStreamingLayer(LSL)
│
├── desktopShortcuts <- Since everyone has to log in to the laptop with their own account, 
│   │                   I gathered all shortcuts for the commonly used software so you can  
│   │                   simply drag these onto your desktop to use them.
│
├── software         <- All the nfb-related software I installed on the laptop (openBCI, 
│   │                   LSL, Arduino, VSCode, miniconda, git, drivers, ...)
│
├── src              <- Source code for this project
│   ├── arduino      <- Script 'nanoWifiMotor.ino' to setup the Arduino Nano as a server to 
│   │                   communicate with a Wifi hotspot on the laptop to use the wristband 
│   │                   for haptic feedback
│   ├── conda        <- Conda environment 'nfb.yml' in which all the python scripts can be run
│   ├── python       <- Python scripts for beta ERD simulation and nfb real-time processing
│   │   │── sim_betaERD   <- simulate beta ERD embedded in pink noise in one channel and 
│   │   │               stream 16 channels to LSL
│   │   │── nfb_fir  <- neurofeedback real-time processing using a fir filter. I tested a lot
│   │   │               but only included this script here to avoid confusion. Also, there 
│   │   │               are improvements I'll make when I find the time, feel free to ask.
├── README.md        <- Top-level README
```