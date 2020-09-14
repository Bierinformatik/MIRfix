import logging
import multiprocessing
import os, sys, inspect
import traceback as tb
import datetime
import shutil
#from MIRfix import scriptname


scriptname = 'MIRfix'


def makelogdir(logdir):
    if not os.path.isabs(logdir):
        logdir =  os.path.abspath(logdir)
    if not os.path.exists(logdir):
        try:
            os.makedirs(logdir)
        except OSError:
            # If directory has already been created or is inaccessible
            if not os.path.exists(logdir):
                sys.exit('Problem creating directory '+logdir)
            else:
                return logdir
    return logdir

def setup_logger(name, log_file, filemode='a', logformat=None, datefmt=None, level='WARNING', proc=1):
    """Function setup as many loggers as you want"""

    if proc > 1:
        log = multiprocessing.get_logger()  # does not take name argument
    else:
        log = logging.getLogger(name)
    if log_file != 'stderr':
        handler = logging.FileHandler(log_file, mode=filemode)
    else:
        handler = logging.StreamHandler(sys.stderr)

    handler.setFormatter(logging.Formatter(fmt=logformat,datefmt=datefmt))

    log.setLevel(level)
    log.addHandler(handler)

    return log

def setup_multiprocess_logger(log_file, filemode='a', logformat=None, datefmt=None, level='WARNING'):
    """Function setup as many loggers as you want"""

    log = multiprocessing.get_logger() # does not take name argument

    if log_file != 'stderr':
        handler = logging.FileHandler(log_file, mode=filemode)
    else:
        handler = logging.StreamHandler(sys.stderr)

    handler.setFormatter(logging.Formatter(fmt=logformat,datefmt=datefmt))

    log.setLevel(level)
    log.addHandler(handler)

    return log

def checklog():
    test = logging.getLogger()
    if not (test.hasHandlers()):
        return False
    else:
        if not len(test.handlers) > 1:
            return False
        else:
            return True

def backup(file):
    if os.path.exists(file):
        os.rename(file,file+'.bak')
    logdir =  os.path.abspath('LOGS')
    if not os.path.exists(logdir):
        os.makedirs(logdir)
        open(os.path.abspath(file),'a').close()

# SETUP log
logtime = str(datetime.datetime.now().strftime("%Y%m%d_%H_%M_%S_%f"))
logdir = 'LOGS' # +os.sep+logtime
makelogdir(logdir)
logfile = os.path.abspath(logdir)+os.sep+scriptname+'_'+logtime+'.log'
log = multiprocessing.get_logger()  # does not take name argument

if scriptname:
    if not os.path.isfile(os.path.abspath(logfile)):
        open(logfile,'a').close()
    else:
        ts = str(datetime.datetime.fromtimestamp(os.path.getmtime(os.path.abspath(logfile))).strftime("%Y%m%d_%H_%M_%S"))
        shutil.move(logfile,logdir+os.sep+scriptname+'_'+ts+'.log')

    setup_multiprocess_logger(log_file='stderr', logformat='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M')
    setup_multiprocess_logger(log_file=logfile, filemode='a', logformat='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M')


if __name__ == '__main__':
    try:
        # set up logging to file
        log = setup_logger(name='', log_file='stderr', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='WARNING')

        # define a Handler which writes INFO messages or higher to the sys.stderr
        #console = logging.StreamHandler()
        #console.setLevel(logging.INFO)
        # set a format which is simpler for console use
        #formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
        # tell the handler to use this format
        #console.setFormatter(formatter)
        # add the handler to the root logger
        #logging.getLogger('').addHandler(console)

        # Now, we can log to the root logger, or any other logger. First the root...
        #logging.info('Imported logger.py')
        # Now, use this in code defining a couple of other loggers which might represent areas in your
        # application, e.g.:
        #log = logging.getLogger('logger.main')

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        logging.error(''.join(tbe.format()))


#def eprint(log, *args, **kwargs):
#    log.error(*args, **kwargs)
#
# log.py ends here
