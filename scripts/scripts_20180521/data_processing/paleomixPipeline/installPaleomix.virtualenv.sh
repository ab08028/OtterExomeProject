#INSTALL PALEOMIX Version: 1.2.13.2
# first install virtualenv
# https://github.com/MikkelSchubert/paleomix/blob/master/docs/installation.rst
pip -t ./ install virtualenv # installed it into my bin dir /u/home/a/ab08028/klohmueldata/annabel_data/bin
python virtualenv.py ~/install/virtualenvs/paleomix
# activate virtual environment
source ~/install/virtualenvs/paleomix/bin/activate
# Following succesful completion of these commands, the paleomix tools will be accessible in the ~/install/virtualenvs/paleomix/bin/ folder. However, as this folder also contains a copy of Python itself, it is not recommended to add it to your PATH. Instead, simply link the paleomix commands to a folder in your PATH. This can, for example, be accomplished as follows:
pip install paleomix ## success!
deactivate
binDir=/u/home/a/ab08028/klohmueldata/annabel_data/bin/paleomixLinks
mkdir -p $binDir
ln -s ~/install/virtualenvs/paleomix/bin/paleomix $binDir
ln -s ~/install/virtualenvs/paleomix/bin/bam_pipeline $binDir
ln -s ~/install/virtualenvs/paleomix/bin/conv_gtf_to_bed $binDir
ln -s ~/install/virtualenvs/paleomix/bin/phylo_pipeline $binDir
ln -s ~/install/virtualenvs/paleomix/bin/bam_rmdup_collapsed $binDir
ln -s ~/install/virtualenvs/paleomix/bin/trim_pipeline $binDir
# it worked! 
# add to path (don't add the virtualenv, add the link dir): export PATH=$PATH:/u/home/a/ab08028/klohmueldata/annabel_data/bin/paleomixLinks
# # Want to test if paleomix works with the shell. Trying making a test file to see
# The test failed. Couldn't find python
# Adding:
# source /u/local/Modules/default/init/modules.sh
# module load python/2.7
# to the script
# THIS WORKS! Note you need to load python module before any paleomix script.

