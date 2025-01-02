subject = 19

n_trialsBlock = 30
n_trialsTotal = 240 + n_trialsBlock # to account for first training block
n_blocks = int(n_trialsTotal/n_trialsBlock)
bslBlocks = 3 # baseline every how many blocks
iti = 0
bslDur = 30 # duration in s
n_corr_train = 10 # 10 # number of consecutive correct trials required for the first training block
audiFeedback = 0 # 1 yes, 0 no
keyList = ['a', 'w', 'e', 'r']
maxRT = 1
maxMT = 3
waitGreen = 1
feedbackDur = 2
use_photosensor = 1 # add a white corner for trigger synchronisation

# Import -------------------------------------------------------------------------------------
from psychopy import core, event, sound
from psychopy.visual import Window, ShapeStim, Rect, TextStim, ImageStim
from pylsl import StreamInfo, StreamOutlet
import random
import pandas
import numpy as np
import os

clock_main = core.Clock()
start_exp = clock_main.getTime()
event_log = []

# Functions ----------------------------------------------------------------------------------
# Function to send LSL markers and save events
def send_triggers(trigger, duration, exp_onset, rt, trial, corr_trial, block):
    outlet.push_sample([trigger])
    
    global clock_main
    global start_exp
    global event_log

    onset = clock_main.getTime() - start_exp

    event_log.append({
        'onset': onset,
        'duration': duration,
        'event': trigger,
        'corr_trial': corr_trial,
        'exp_onset': exp_onset,
        'trial': trial+1,
        'block': block+1,
        'response_time': rt
    })

# Create events DataFrame and save to .tsv file
def save_events(event_log):
    filepath = f'..\\..\\rec\\simpleRT\\sub-P{"%03d" % subject}\\ses-S001\\eeg'
    filename = f'sub-{"%03d" % subject}_ses-001_task-simpleRT_run-001_events.tsv'
    fullfile = os.path.join(filepath, filename)
    os.makedirs(filepath, exist_ok=True) 

    df = pandas.DataFrame(event_log)
    df.to_csv(fullfile, sep="\t", index=False)

def safe_quit(key):
    if key == 'escape':
        save_events(event_log)
        send_triggers('escape', 'n/a', 'n/a', 'n/a', -1, 'n/a', -1)
        win.close()
        core.quit()

# Setup --------------------------------------------------------------------------------------
# Setup LSL
info = StreamInfo(name="LSL_Markers", type="Markers", channel_count=1,
                  channel_format="string", source_id="lsl_markers")
outlet = StreamOutlet(info)

# Setup PsychoPy Window
win_size = 1920, 1080  # window size
win = Window(size=win_size, fullscr=True, allowGUI = False, color=(-1, -1, -1),
                    units='pix', winType='pyglet', screen=1)

# Create a fixation cross
fixCross = ShapeStim(win, vertices=((0, -20), (0, 20), (0,0), (-20,0), (20, 0)),
                            lineWidth=5, closeShape=False, lineColor=[1, 1, 1])

# Create a green square
green_square = Rect(win, width=250, height=250, fillColor=[0, 1, 0], lineColor=[0, 1, 0])

# Progress bar
progress_bar = Rect(win, width=600, height=50, fillColor=None, lineColor=[1, 1, 1], pos=(0, -200))

# Load movement sequence to display at top of screen for learning
scale_img = .5
img_seq = ImageStim(win, image='../../img/nfb_fingerSeq_left.png', pos=(0, 310), size=(832*scale_img, 209*scale_img)) # 832x209 is the image size in pixels

# Small white rectangle to place photosensor on for trigger sync
stim_size = 100, 100  # stimulus size
LOWER_RIGHT = (win_size[0]/2 - stim_size[0]/2, stim_size[1]/2 - win_size[1]/2)
BLACK, WHITE = (-1, -1, -1), (1, 1, 1) # stimulus colors
if use_photosensor:
    photoRectOn = Rect(win=win, width=stim_size[0], height=stim_size[1], fillColor=WHITE, pos=LOWER_RIGHT)
else:
    photoRectOn = Rect(win=win, width=stim_size[0], height=stim_size[1], fillColor=BLACK, pos=LOWER_RIGHT) # stays black throughout
photoRectOff = Rect(win=win, width=stim_size[0], height=stim_size[1], fillColor=BLACK, pos=LOWER_RIGHT)
# Wait for 1 second
core.wait(1)

# Wait for key press to move on so LSL streams can be checked in the recorder
TextStim(win, text='Experimenter:\n\nSet correct subject in VSCode.\nStart LSL Recording.\n\n\n\n\nPress \'space\' to continue', color=(1, 1, 1)).draw()
win.flip()
event.waitKeys(keyList=['space'])

# Count points
points_counter = 0

# First trigger
send_triggers('expStart', 'n/a', 'n/a', 'n/a', -1, 'n/a', -1)

# Main loop around blocks -------------------------------------------------------------------------
for block in range(n_blocks):

    # Comment in to skip blocks that were already done if we had to take a break during the experiment
    # if block == 0 or block == 1: # example to skip the first two blocks
    #     continue

    # Allow for more trials in the first block so the participants can learn the sequence
    if block == 0:
        n_trials = 200
    else:
        n_trials = n_trialsBlock

    # Randomization -------------------------------------------------------------------------------
    rseed = float(f"{subject}.{block}")
    random.seed(subject)

    # Generate 9 equidistant values between 1 and 5
    iti_values = np.linspace(1, 5, 9)

    # Repeat the values enough times
    repeated_iti_values = np.tile(iti_values, (n_trials // len(iti_values)) + 1)

    # Trim the array to the right number of elements
    goCue = repeated_iti_values[:n_trials]

    # Shuffle the values
    np.random.shuffle(goCue)

    # Other option would have been to generate non-equidistant random uniformly distributed iti values:
    # goCue = [round(random.uniform(iti,5),3) for _ in range(n_trials)]

    # Timing --------------------------------------------------------------------------------------
    clock_bsl = core.Clock()
    start_bsl = clock_bsl.getTime()

    # Baseline ------------------------------------------------------------------------------------
    if block == 1:
        eegSetupInstructions =  TextStim(win, text='Setup EEG\n\n\nPress \'space\' to continue', color=(1, 1, 1))
        eegSetupInstructions.draw()
        win.flip()
        event.waitKeys(keyList=['space']) # Wait for the participant to press 'space' to continue
    if block % bslBlocks == 0 or block == 1:
        bslInstructions =  TextStim(win, text='For the next 30 seconds:\n\nMaintain fixation on the cross in the screen center.\n\nPlease rest and stay still!\n\n\n\n\nPress \'space\' to continue', color=(1, 1, 1))
        bslInstructions.draw()
        win.flip()
        event.waitKeys(keyList=['space']) # Wait for the participant to press 'space' to continue
        win.callOnFlip(send_triggers, 'bslStart', bslDur, 'n/a', 'n/a', -1, 'n/a', block)
        fixCross.draw()
        win.flip()
        core.wait(bslDur)
        win.callOnFlip(send_triggers, 'bslStop', 'n/a', 'n/a', 'n/a', -1, 'n/a', block)

    # Experiment instructions ---------------------------------------------------------------------
    if block == 0:
        win.callOnFlip(send_triggers, 'instructions_00', 'n/a', 'n/a', 'n/a', -1, 'n/a', block)
        TextStim(win, text=f'Part 1\n\n\nPlace your left hand on the laptop keyboard like this:\nr: index\ne: middle\nw: ring\na:pinky\n\n\n\nPress \'space\' to continue', color=(1, 1, 1)).draw()
        win.flip()
        event.waitKeys(keyList=['space']) # Wait for the participant to press 'space' to continue

        win.callOnFlip(send_triggers, 'instructions_01', 'n/a', 'n/a', 'n/a', -1, 'n/a', block)
        expInstructions = TextStim(win, text='Part 1\n\n\nMaintain fixation on the cross in the screen center.\n\nWhen you see a green rectangle, carry out this finger sequence: index - ring - middle - pinky\n\n\n\n\n\n\n\n\n\n\n\nTraining Block: The sequence will be displayed at the top of the screen.\nTry to memorize it!\n\nRespond as quickly and accurately as possible!\n\n\n\n\nPress \'space\' to continue', color=(1, 1, 1))

        # Images
        images = ['fixCross.png','goCue.png','index_left.png', 'ring_left.png', 'middle_left.png', 'pinkie_left.png']
        image_paths = ['../../img/digits/' + image for image in images]# Create ImageStim objects
        img_Instr = [ImageStim(win, image=image_path, pos=(0,30)) for image_path in image_paths]
        img_index = 0

        clock_instr = core.Clock()
        while True:
            # Draw the current image and instructions
            img_Instr[img_index].draw()
            expInstructions.draw()
            win.flip()
            
            # Check if 'space' is pressed
            keys = event.getKeys(keyList=['space'])
            if 'space' in keys:
                break
            
            # Cycle through images every 1 s
            if clock_instr.getTime() > 1:
                img_index = (img_index + 1) % len(images)
                clock_instr.reset()
    else:
        win.callOnFlip(send_triggers, f'instructions_{block}0', 'n/a', 'n/a', 'n/a', -1, 'n/a', block)
        TextStim(win, text=f'Part 1\n\n\nCompleted block {block-1} of {n_blocks-1}.\n\n\n\nTake a break for as long as you\'d like.\n\n\n\n\nPress \'space\' when you\'re ready to continue.', color=(1, 1, 1)).draw()
        win.flip()
        event.waitKeys(keyList=['space']) # Wait for the participant to press 'space' to continue
        win.callOnFlip(send_triggers, f'instructions_{block}1', 'n/a', 'n/a', 'n/a', -1, 'n/a', block)
        TextStim(win, text=f'Part 1\n\n\nPlace your left hand on the laptop keyboard like this:\nr: index\ne: middle\nw: ring\na:pinky\n\n\n\nPress \'space\' to continue.', color=(1, 1, 1)).draw()
        win.flip()
        event.waitKeys(keyList=['space'])
        win.callOnFlip(send_triggers, f'instructions_{block}2', 'n/a', 'n/a', 'n/a', -1, 'n/a', block)
        TextStim(win, text=f'Part 1\n\n\nBlock {block} / {n_blocks-1}\n\n\nMaintain fixation on the cross in the screen center.\n\n\nRespond as quickly and accurately as possible!\n\n\n\n\nPress \'space\' to start the block.', color=(1, 1, 1)).draw()
        win.flip()
        event.waitKeys(keyList=['space'])

    # Main experiment loop -------------------------------------------------------------------------    
    clock_block = core.Clock()

    fixCross.draw()
    start_time = clock_block.getTime()
    send_triggers('blockStart', goCue[0], 0, 'n/a', 0, 'n/a', block)
    win.flip()

    prev_start_time = start_time
    corr_counter = 0

    for trial in range(n_trials):

        # Calculate flip_time for this trial
        flip_time = prev_start_time + goCue[trial]

        # Display sequence at top of screen for learning during first block
        if block == 0  and corr_counter < n_corr_train/2: 
            img_seq.draw()

        # Draw fixation cross and green square
        green_square.draw()
        photoRectOn.draw()
        fixCross.draw()

        # Ensure proper timing of stimulus onsets
        while clock_block.getTime() < (flip_time - .003): # allow it to move on x ms before flip time so it catches the next flip. Otherwise there's a constant offset of ~14 ms.
            keys = event.getKeys(timeStamped=clock_block)
            for key, key_time in keys:
                send_triggers(key, 1, 'n/a', key_time - (prev_start_time - 3 - iti), trial-1, 'n/a', block)
                safe_quit(key)
        win.flip()
        send_triggers('go', maxRT + waitGreen, flip_time, 'n/a', trial, 'n/a', block)

        # Listen for key presses
        key_presses = []
        rt_sequence = []
        first_rt = None
        while clock_block.getTime() - flip_time < maxMT + waitGreen - .3: 
            keys = event.getKeys(timeStamped=clock_block)
            for key, key_time in keys:
                send_triggers(key, 1, 'n/a', key_time - flip_time, trial, 'n/a', block)
                if (key_time  - flip_time) < .02:
                    continue
                if first_rt is None:
                    first_rt = key_time - flip_time
                key_presses.append(key)
                rt_sequence.append(key_time - flip_time)
                safe_quit(key)
        
        # Check correctness
        correct_sequence = ['r', 'w', 'e', 'a']

        if len(key_presses) < 4 or len(key_presses) > 4 or key_presses[:4] != correct_sequence:
            feedback = 'Error: Incorrect sequence'
            corr_trial = 0
            corr_counter = 0
            feedback_tone = sound.Sound('C', secs=0.1)  # Tone C (261.63 Hz)
        elif first_rt is not None and first_rt > maxRT:
            feedback = 'Error: Initial reaction time too slow'
            corr_trial = 0
            corr_counter = 0
            feedback_tone = sound.Sound('D', secs=0.1) 
        elif rt_sequence[3] > maxMT:
            feedback = 'Error: Movement sequence too slow'
            corr_trial = 0
            corr_counter = 0
            feedback_tone = sound.Sound('E', secs=0.1)
        else:
            corr_trial = 1
            corr_counter += 1
            points_counter += 10
            if block == 0:
                feedback = f'Correct trial'
            else:
                feedback = f'\n\n\n\nCorrect trial\n\n\n\n+10 points = {points_counter} points'
            feedback_tone = sound.Sound('A', secs=0.1)  # Tone A (440 Hz)

        # Display feedback
        # Progress bar
        filled_width = (corr_counter / n_corr_train) * 600  # 600 is the full width of the progress bar
        progress_bar_filled = Rect(win, width=filled_width, height=50, fillColor=[0, 1, 0], lineColor=[0, 1, 0], pos=(filled_width/2 - 300, -200))
        # Message
        feedback_message = TextStim(win, text=feedback, color=(1, 1, 1))
        feedback_message.draw()
        if block == 0: 
            progress_bar.draw()  # Draw the full progress bar outline
            progress_bar_filled.draw()  # Draw the filled part of the progress bar
            if corr_counter < n_corr_train/2:
                img_seq.draw()
        while clock_block.getTime() < (flip_time + maxMT + waitGreen - .002): # leave 2ms to spare again
            pass
        photoRectOff.draw()
        win.flip()
        send_triggers('feedback', feedbackDur, flip_time + maxMT + waitGreen, 'n/a', trial, corr_trial, block)

        # Auditory feedback
        if audiFeedback:
            feedback_tone.play()
            core.wait(.1)
            feedback_tone.stop() # Otherwise the feedback tone can't be played again

        # Quit if we have n_corr_train consecutive correct trials in the training block
        if block == 0 and corr_counter == n_corr_train:
            points_counter = 0 # reset after training for first experimental block
            core.wait(feedbackDur)
            break

        # Display only fixation cross during iti
        fixCross.draw()
        if block == 0 and corr_counter < n_corr_train/2: 
            img_seq.draw() # Display sequence at top of screen for learning during first block
        while clock_block.getTime() < (flip_time + maxMT + waitGreen + feedbackDur - .002): # leave 2 ms to spare again
            keys = event.getKeys(timeStamped=clock_block)
            for key, key_time in keys:
                send_triggers(key, 1, 'n/a', key_time - flip_time, trial, 'n/a', block)
                safe_quit(key)
        win.flip()
        if trial < n_trials - 1:
            send_triggers('iti', iti + goCue[trial + 1], flip_time + maxMT + waitGreen + feedbackDur, 'n/a', trial, 'n/a', block)
        else:
            send_triggers('iti', 0, flip_time + maxMT + waitGreen + feedbackDur, 'n/a', trial, 'n/a', block)

        # Update start_time for the next trial
        prev_start_time = flip_time + maxMT + waitGreen + feedbackDur + iti

        # Save events
        save_events(event_log)

# Baseline period
bslInstructions =  TextStim(win, text='For the next 30 seconds:\n\nMaintain fixation on the cross in the screen center.\n\nRelax your mind and don\'t move!\n\n\n\n\nPress \'space\' to continue', color=(1, 1, 1))
bslInstructions.draw()
win.flip()
event.waitKeys(keyList=['space'])
win.callOnFlip(send_triggers, 'bslStart', bslDur, 0, 'n/a', -1, 'n/a', -1)
fixCross.draw()
win.flip()
core.wait(bslDur)
send_triggers('bslStop', 'n/a', 0, 'n/a', -1, 'n/a', -1)

# End of experiment
bslInstructions =  TextStim(win, text='Experiment finished\n\nWell done!!!\n\n\n\n\nPress \'space\' to close.', color=(1, 1, 1))
bslInstructions.draw()
win.flip()
send_triggers('expStop', 'n/a', 'n/a', 'n/a', -1, 'n/a', -1)
event.waitKeys(keyList=['space']) # Wait for the participant to press "p" to continue

# Save .tsv file with events
save_events(event_log)

win.close()
core.quit()
