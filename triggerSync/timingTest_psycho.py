# Import
from psychopy import core, event
from psychopy.visual import Window, Rect
from pylsl import StreamInfo, StreamOutlet
import serial
import random as rand
#import pyglet

def send_triggers(value):
    """Send value as hardware trigger and LSL marker."""
    outlet.push_sample(value)
    # port.write('WRITE {}\n'.format(str(value)[1]).encode())

# Setup LSL
# Create stream info
info = StreamInfo(name="LSL_Markers", type="Markers", channel_count=1,
                  channel_format="int32", source_id="LSL_Markers_00{}".format(rand.randint(1,9)))
# Make an outlet
outlet = StreamOutlet(info)

# Setup ParallelPort
# port = serial.Serial("/dev/ttyACM0", baudrate=128000, bytesize=8)

# Settings

BLACK, WHITE = (-1, -1, -1), (1, 1, 1) # stimulus colors

n_trials = 1000 # number of trials
win_size = 1080, 1920  # window size
stim_size = 100, 100  # stimulus size
n_frames = 1  # number of frames: 1 frame at 60 Hz refresh rate is 16 ms

win = Window(size=win_size, fullscr=True, allowGUI = False, color=BLACK,
             units='pix', winType='pyglet', screen=1)

# Stimulus position
LOWER_LEFT = (stim_size[0]/2 - win_size[0]/2, stim_size[1]/2 - win_size[1]/2)
UPPER_LEFT = (stim_size[0]/2 - win_size[0]/2, win_size[1]/2 - stim_size[1]/2)
LOWER_RIGHT = (win_size[0]/2 - stim_size[0]/2, stim_size[1]/2 - win_size[1]/2)
UPPER_RIGHT = (win_size[0]/2 - stim_size[0]/2, win_size[1]/2 - stim_size[1]/2)

rect1 = Rect(win=win, width=stim_size[0], height=stim_size[1], fillColor=WHITE, pos=LOWER_RIGHT)
rect2 = Rect(win=win, width=stim_size[0], height=stim_size[1], fillColor=BLACK, pos=LOWER_RIGHT)

# Monitor refresh rate
refreshRate = win.getMsPerFrame(nFrames=240) # Get monitor refresh rate

# Stop script if monitor refresh rate is not set to 60 Hz
if refreshRate[2] <  16 or refreshRate[2] >  17:
    print("Refresh rate has to be 60 Hz!")
    win.close()
    core.quit()

# Main loop ----------------------------------------------------------
core.wait(2) # Wait before paradigm starts 

# Send triggers every time an image changes
for i in range(n_trials):
    if event.getKeys(keyList=["escape"]):
        break
    win.callOnFlip(send_triggers, [1])
    for n in range(n_frames):
        rect1.draw()
        win.flip()
    win.callOnFlip(send_triggers, [2])
    for n in range(n_frames):
        rect2.draw()
        win.flip()
    core.wait(0.1+rand.random()/1000*40) # jitter onsets uniform 80-120 ms 

# Close & quit
# port.close()
win.close()
core.quit()