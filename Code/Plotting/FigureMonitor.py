import time
from watchdog.observers import Observer
from watchdog.events import PatternMatchingEventHandler
import sys
from PlotFrameFunction import graph
import os
#shamelessly from http://thepythoncorner.com/dev/how-to-create-a-watchdog-in-python-to-look-for-filesystem-changes/

cmdargs = list(map(str,sys.argv))

mdirectory =  cmdargs[1]

colarray='red'


if(len(cmdargs)> 2):
    x= False
    colarray = cmdargs[2]
else:
    x= True
    

n=0

def myfunc() :
    global n
    n = 0

if __name__ == "__main__":
    patterns = "*"
    ignore_patterns = ""
    ignore_directories = True
    case_sensitive = True
    my_event_handler = PatternMatchingEventHandler(patterns, ignore_patterns, ignore_directories, case_sensitive)

def on_created(event):
    #print(f"hey, {event.src_path} has been created!")
    myfunc()
    filecreated = str(event.src_path)
    filename, file_extension = os.path.splitext((filecreated))
    typ = filename[2:5]
    if file_extension == ".csv" and typ == "pos":
        graph(filecreated,colarray)
    else:
        print("not csv created")
    # if x>0:
    #     
    # else:
    #     print("not csv created")

    #if({event.src_path}.endswith)
    #

def on_deleted(event):
    myfunc()
#    print(f"Someone deleted {event.src_path}!")

def on_modified(event):
    myfunc()
#    print(f"hey buddy, {event.src_path} has been modified")

def on_moved(event):
    myfunc()
#    print(f"ok ok ok, someone moved {event.src_path} to {event.dest_path}")

my_event_handler.on_created = on_created
my_event_handler.on_deleted = on_deleted
my_event_handler.on_modified = on_modified
my_event_handler.on_moved = on_moved


path = mdirectory
go_recursively = True
my_observer = Observer()
my_observer.schedule(my_event_handler, path, recursive=go_recursively)

my_observer.start()

try:
    while (n<100):
        time.sleep(1)
        n=n+1
        #print(n)
except KeyboardInterrupt:
    my_observer.stop()
    my_observer.join()    