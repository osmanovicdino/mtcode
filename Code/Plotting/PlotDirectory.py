import time
from watchdog.observers import Observer
from watchdog.events import PatternMatchingEventHandler
import sys
from PlotFrameFunction import graph
import os
#shamelessly from http://thepythoncorner.com/dev/how-to-create-a-watchdog-in-python-to-look-for-filesystem-changes/

cmdargs = list(map(str,sys.argv))

mdirectory =  cmdargs[1]


colarraytemp='red'

