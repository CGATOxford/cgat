import glob

__all__ = [x[:-3] for x in glob.glob("*.py")]
