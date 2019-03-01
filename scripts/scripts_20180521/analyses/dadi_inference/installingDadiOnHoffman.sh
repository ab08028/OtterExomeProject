#installing dadi on hoffman: 
#https://support.idre.ucla.edu/helpdesk/KB/View/11345064-python-environments-on-hoffman
qrsh 
module load python/2.7.13_shared
virtualenv $HOME/env_python2.7.13
source $HOME/env_python2.7.13/bin/activate
pip install numpy
pip install scipy
cd /u/home/a/ab08028/bin
git clone https://bitbucket.org/gutenkunstlab/dadi
cd dadi 
# build dadi 1.6.3
python setup.py install # installs it just in this virutalenv
# YOU MUST LEAVE THE DADI DIRECTORY TO IMPORT DADI! <-- OTHERWISE YOU GET THE TRIDIAG ERROR
cd ../
# test import of dadi in python ; python --> import dadi  # (works)
# also need matplotlib
deactivate # deactivates your virtual env


# Then when you use dadi in a job you need to do this at the start of your bash submission script: 
source /u/local/Modules/default/init/modules.sh
module load python/2.7.13_shared
source /u/home/a/ab08028/env_python2.7.13/bin/activate # activate virtual environment for dadi

# then at the end of the script:
deactivate