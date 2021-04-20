# delete stuff in build dir
cd ../cvode-5.7.0
rm -r build_dir

# yaml dir
cd ../fortran-yaml
rm -r build_dir

cd ../
chmod -R +w StringiFor
rm -r StringiFor

cd dependencies
find . -type f ! -name '*.sh' -delete
rm -r -- ./*/
