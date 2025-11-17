import time
from contextlib import contextmanager

# --- Setting up timing
@contextmanager
def log_timed(logger, msg):
    """Time a piece of code in a with context, yielding a message with time spent in seconds.
    
    Parameters
    ----------
    logger: logging.logger
        (local logger obtained from logging.getlogger())
    msg: str
        the message.
        
    Usage
    -----
    with log_time(logger, msg):
        ...
        ...
    
    """
    start = time.perf_counter()
    yield
    logger.info(f"{msg} in {time.perf_counter() - start:.2f} seconds")
