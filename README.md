23.02.2014

Should really fix the fact that the model radius array always needs to have a minus 1 after the length.
Shouldn't be too hard to do.

11.02.2014

The function 'to_string' for some reason is not supported with the intel standard libraries. To get this working on monch you need to
	'module load gcc'
in addition to the intel libraries. Pretty dumb.
