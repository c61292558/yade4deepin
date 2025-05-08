"Module providing execfile forward compatibility with new python versions"


def execfile(filename, myglobals=None, mylocals=None):
	"""
        Read and execute a Python script from a file in the given namespaces.
        The globals and locals are dictionaries, defaulting to the current
        globals and locals. If only globals is given, locals defaults to it.
        """
	try:
		from collections.abc import Mapping
	except:
		from collections import Mapping
	import inspect
	if myglobals is None:
		# There seems to be no alternative to frame hacking here.
		caller_frame = inspect.stack()[1]
		myglobals = caller_frame[0].f_globals
		mylocals = caller_frame[0].f_locals
	elif mylocals is None:
		# Only if myglobals is given do we set mylocals to it.
		mylocals = myglobals
	if not isinstance(myglobals, Mapping):
		raise TypeError('globals must be a mapping')
	if not isinstance(mylocals, Mapping):
		raise TypeError('locals must be a mapping')
	with open(filename, "rb") as fin:
		source = fin.read()
	code = compile(source, filename, "exec")
	## TODO: maybe try deepcopy of myglobals, mylocals so that execfile does not modify the original globals and locals
	exec(code, myglobals, mylocals)
