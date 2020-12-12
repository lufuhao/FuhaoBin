#!/usr/bin/env python3

###For pip >= 10.0.1
import pkg_resources
from subprocess import call

packages = [dist.project_name for dist in pkg_resources.working_set]
call("pip install --upgrade " + ' '.join(packages), shell=True)


'''
### Method2
pip list --outdated --format=freeze | grep -v '^\-e' | cut -d = -f 1  | xargs -n1 pip install -U


### pip < 10.0.1
import pip
from subprocess import call
 
for dist in pip.get_installed_distributions():
    call("pip install --upgrade " + dist.project_name, shell=True)
'''
