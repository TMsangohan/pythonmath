#!usr/bin/python

# Example to run bash script from python

import subprocess

subprocess.call("script.sh",shell=True)

# when using providing input to the bash script

Process=Popen('./childdir/execute.sh %s %s' % (str(var1),str(var2),), shell=True)