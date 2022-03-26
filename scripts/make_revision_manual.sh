######################################################################
# unpacks documentation tarfile to the homepage
#
#    bash make_revision_manual.sh 6 3 0
#
# to generate documentation for release 6.3.0
######################################################################

MAJOR=$1
MINOR=$2
PATCH=$3
DIR=$MAJOR.$MINOR.$PATCH

cd /var/www/html
mkdir $DIR
cd $DIR
tar xvzf ../download/deal.doc-$DIR.tar.gz
mv deal.II/doc/* .
rmdir deal.II/doc
rmdir deal.II
chmod -R a+r .
for i in `find .` ; do if test -d $i ; then chmod a+x $i ; fi ; done
