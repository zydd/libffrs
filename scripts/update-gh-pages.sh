set -euxo pipefail

git checkout master
git pull

DESCRIPTION=$(git describe)

source .venv/bin/activate

rm -rf build/html/
rm -rf docs/generated/

make -C docs/ html

git checkout gh-pages
git pull

rm -rf *.html *.inv *.js
rm -rf _sources/ extra/ _static/ generated/

cp -r build/html/* .

git add *.html *.inv *.js
git add _sources/ extra/ _static/ generated/

git commit -a -m"Update to $DESCRIPTION"
