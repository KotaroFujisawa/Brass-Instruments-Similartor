"""Render the double horn's four-octave concert B-flat major scale."""
import sys

from render_euphonium_ophicleide import main


if __name__=='__main__':
    main(['horn',*sys.argv[1:]])
