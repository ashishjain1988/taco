wget https://sourceforge.net/projects/rna-cpat/files/v1.2.2/CPAT-1.2.2.tar.gz
tar -zxvf CPAT-1.2.2.tar.gz
cd CPAT-1.2.2
mkdir -p cpat_build/lib/python2.7/site-packages
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CPATBUILDDIR=$DIR/cpat_build
export PYTHONPATH=$CPATBUILDDIR/lib/python2.7/site-packages
sed -e '20,21d' setup.py > setup.clip.py
python setup.clip.py install --prefix=$CPATBUILDDIR
pyinstaller --onefile -n cpat_exec --hidden-import numpy --hidden-import TabProxies  build/scripts-2.7/cpat.py
