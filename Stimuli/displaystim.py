import StimulusRoutines as stim
from DisplayStimulus import DisplaySequence
flashing_circle = stim.FlashingCircle(mon,ind)
ds = DisplaySequence(log_dir='log_directory')
ds.trigger_display()
