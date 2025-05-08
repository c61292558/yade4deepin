# jerome.duriez@inrae.fr
def ramFP():
	'''Returns the RAM footprint in MB = 1024 kB'''
	import resource
	mem = resource.getrusage(resource.RUSAGE_SELF)[
	        2
	]  # maximum resident set size used in kB (see ru_maxrss https://manpages.debian.org/buster/manpages-dev/getrusage.2.en.html) = RAM footprint (without swap) https://en.wikipedia.org/wiki/Resident_set_size
	return mem / 1024.  # returns MB
