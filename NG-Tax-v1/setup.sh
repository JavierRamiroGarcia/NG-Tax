sudo apt-get install dos2unix libfreetype6-dev python-pip git  python-numpy
git clone https://github.com/biocore/qiime.git
cd qiime
#stick to certain checkout so we are sure its keeps working and the result do not change overtime
git checkout 1.9.0
pip install .
echo "y" | pip uninstall matplotlib
pip install matplotlib==1.2.0
cd ..
echo -e "PATH=~/.local/bin/:$PATH" >> ~/.bashrc
source ~/.bashrc


