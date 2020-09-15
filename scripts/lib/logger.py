import logging
import multiprocessing
import os
import sys
import inspect
import traceback as tb
import datetime
import shutil


scriptname = __name__


class logger():
    def __init__(self, name, logdir=None, logfile=None, level=None, lfmt='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', dfmt='%m-%d %H:%M', multi=None):
        if multi:
            self.logger = multiprocessing.Logger()
        else:
            self.logger = logging.getLogger(name)

        self.logtime = str(datetime.datetime.now().strftime("%Y%m%d_%H_%M_%S_%f"))
        self.logdir = logdir
        self.logfile = logfile
        self.logformat = lfmt
        self.datefmt = dfmt
        self.propagate = True
        self.stream_handlers = None
        self.logfile_handlers = None

    def _makelogdir(self):
        if not os.path.isabs(self.logdir):
            self.logdir = os.path.abspath(self.logdir)
        if not os.path.exists(self.logdir):
            try:
                os.makedirs(self.logdir)
            except OSError:
                # If directory has already been created or is inaccessible
                if not os.path.exists(self.logdir):
                    sys.exit('Problem creating directory '+self.logdir)


    def _setup_logger(self):
        """Function setup as many loggers as you want"""
        self.stream_handler = logging.StreamHandler(sys.stderr)
        self.stream_handler.setFormatter(logging.Formatter(fmt=self.logformat, datefmt=self.datefmt))
        self.logger.addHandler(self.stream_handler)


    def _addHandler(self):
        if self.logfile != 'stderr':
            self.logfile_handler = logging.FileHandler(self.logfile, mode='a')
            self.logfile_handler.setFormatter(logging.Formatter(fmt=self.lfmt, datefmt=self.dfmt))
            self.logger.addHandler(self.logifle_handler)


    def _checklog(self):
        test = logging.getLogger()
        if not (test.hasHandlers()):
            return False
        else:
            if not len(test.handlers) > 1:
                return False
            else:
                return True


    def _makelogfile(self):
        if not os.path.isfile(os.path.abspath(self.logfile)) or os.stat(self.logfile).st_size == 0:
            open(self.logfile,'a').close()
        else:
            ts = str(datetime.datetime.fromtimestamp(os.path.getmtime(os.path.abspath(self.logfile))).strftime("%Y%m%d_%H_%M_%S"))
            shutil.move(self.logfile,logfile.replace('.log', '')+'_'+ts+'.log')


    def info(self, msg):
        self.logger.info(msg)


    def debug(self, msg):
        self.logger.debug(msg)


    def warning(self, msg):
        self.logger.warning(msg)


    def error(self, msg):
        self.logger.error(msg)


if __name__ == '__main__':
    try:
        # set up logging to file
        logger = logger()
        logger.setup_logger(name='', log_file='stderr', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='WARNING')

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
